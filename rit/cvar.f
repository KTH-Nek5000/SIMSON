c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cvar(ur,ivar,pxz,pxz2,pxz3,pxy,pxy2,w,wxz,
     &   prey,prex,prez,re,eta,alfa,beta,dstar,fltype,
     &   phiuw,uwsub)

      implicit none

      include 'par.f'

      complex ur(memnx,memny,memnz,5)
      complex pxy(nx/2,nyp),pxy2(nx/2,nyp),w(nx/2,nyp)
      complex pxy3(nx/2,nyp)
      real pxz(nx+2,nz),pxz2(nx+2,nz),pxz3(nx+2,nz),wxz(nx+2,nz)
      real prey(nyp*2+15),alfa(nx/2),beta(nz),prex(nx+15),prez(nz*2+15)
      real re,dstar,phiuw
      integer ivar
      logical uwsub
      integer fltype
      real eta(nyp)

      integer x,y,z,i
      real sym,tr
c
c     omx=wy-vz
c
      if (ivar.eq.4) then
         do z=1,nzc
            call getxyp(pxy,z,3,ur)
            call vchbf(pxy,w,nyp,nx,nx,1,prey)
            call rdcheb(pxy,nyp,nx,nx)
            if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
               do y=1,nyp
                  do x=1,nx/2
                     pxy(x,y)=pxy(x,y)*dstar
                  end do
               end do
            end if
            call vchbb(pxy,w,nyp,nx,nx,1,prey)
            call getxyp(pxy2,z,2,ur)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=pxy(x,y)*(2./real(nyp-1))-(0.,1.)*beta(z)
     &                 *pxy2(x,y)
               end do
            end do
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     omy=uz-wx
c
      if (ivar.eq.5) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            call getxyp(pxy2,z,3,ur)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=(0.,1.)*beta(z)*pxy(x,y)-(0.,1.)*alfa(x)*
     &                 pxy2(x,y)
               end do
            end do
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     omz=vx-uy
c
      if (ivar.eq.6) then
         do z=1,nzc
            call getxyp(pxy,z,2,ur)
            call getxyp(pxy2,z,1,ur)
            call vchbf(pxy2,w,nyp,nx,nx,1,prey)
            call rdcheb(pxy2,nyp,nx,nx)
            if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
               do y=1,nyp
                  do x=1,nx/2
                     pxy2(x,y)=pxy2(x,y)*dstar
                  end do
               end do
            end if
            call vchbb(pxy2,w,nyp,nx,nx,1,prey)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=(0.,1.)*alfa(x)*pxy(x,y)-pxy2(x,y)*
     &                 (2./real(nyp-1))
               end do
            end do
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     u-umean
c
      if (ivar.eq.7) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            if (z.eq.1) then
               do y=1,nyp
                  pxy(1,y)=(0.0,0.0)
               end do
            end if
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     u-ulam
c
      if (ivar.eq.8) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            if (z.eq.1) then
               do y=1,nyp
                  if (fltype.eq.1.or.fltype.eq.4) then ! Needs to be checked???
                     pxy(1,y)=pxy(1,y)-1.+eta(y)**2-pxy(1,nyp)
                  else
                     pxy(1,y)=pxy(1,y)-eta(y)
                  end if
               end do
            end if
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     dudy
c
      if (ivar.eq.9.or.ivar.eq.14) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            call vchbf(pxy,w,nyp,nx,nx,1,prey)
            call rdcheb(pxy,nyp,nx,nx)
            if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
               do y=1,nyp
                  do x=1,nx/2
                     pxy(x,y)=pxy(x,y)*dstar
                  end do
               end do
            end if
            call vchbb(pxy,w,nyp,nx,nx,1,prey)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=pxy(x,y)*(2./real(nyp-1))
               end do
            end do
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     Energy in physical space
c
      if (ivar.eq.10) then
         do y=1,nyp
            call getxzp(pxz,y,1,ur,sym) ! get u
            call vcfftb(pxz,pxz(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz,pxz(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call getxzp(pxz2,y,2,ur,sym) ! get v
            call vcfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call getxzp(pxz3,y,3,ur,sym) ! get w
            call vcfftb(pxz3,pxz3(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz3,pxz3(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)

            do z=1,nzc
               do x=1,nx
                 pxz(x,z)=(pxz(x,z)**2+pxz2(x,z)**2+pxz3(x,z)**2)/
     &                     real(2*nx*nz)
               end do
            end do

            call vrfftf(pxz,pxz(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call vcfftf(pxz,pxz(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call putxzp(pxz,y,4,ur)
         end do
      end if
c
c      Spectral energy, note that a square root is taken to ensure
c      compatibility with routines calculating spectra of a velocity
c      this function was originally in Chebyshev spectral space
c      but this was not useful for certain purposes since
c      the Chebyshev polynomials are not orthogonal with weight 1
c      this was changed 910430.
c      the old function is now available for ivar = 12
c
      if (ivar.eq.11.or.ivar.eq.12) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            if (ivar.eq.12) call vchbf(pxy,w,nyp,nx,nx,1,prey)
            do y=1,nyp
               do x=1,nx/2
                  pxy2(x,y)=pxy(x,y)*conjg(pxy(x,y))
               end do
            end do
            do i=2,3
               call getxyp(pxy,z,i,ur)
               if (ivar.eq.12)  call vchbf(pxy,w,nyp,nx,nx,1,prey)
               do y=1,nyp
                  do x=1,nx/2
                     pxy2(x,y)=pxy2(x,y)+pxy(x,y)*conjg(pxy(x,y))
                  end do
               end do
            end do
            do y=1,nyp
               do x=1,nx/2
                  if (ivar.eq.12) then
                     pxy2(x,y)=sqrt(real(pxy2(x,y)))*2./real(nyp-1)
                  else
                     pxy2(x,y)=sqrt(real(pxy2(x,y)))
                  end if
               end do
            end do
            if (ivar.eq.12) call vchbb(pxy2,w,nyp,nx,nx,1,prey)
            call putxyp(pxy2,z,4,ur)
         end do
      end if
c
c      Reynolds stress = -(u-umean)*v or
c      Total shear = -(u-umean)*v+dudy*nu
c      or turbulence production -(u-umean)*v*dumeandy
c      Note that dudy is already calculated above
c
      if (ivar.ge.13.and.ivar.le.15) then
         if (ivar.eq.15) then
            call getxyp(pxy,1,1,ur)
            call vchbf(pxy,w,nyp,nx,nx,1,prey)
            call rdcheb(pxy,nyp,nx,nx)
            if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
               do y=1,nyp
                  do x=1,nx/2
                     pxy(x,y)=pxy(x,y)*dstar
                  end do
               end do
            end if
            call vchbb(pxy,w,nyp,nx,nx,1,prey)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=pxy(x,y)*(2./real(nyp-1))
               end do
            end do
         end if
         do y=1,nyp
            sym=1.
            call getxzp(pxz,y,1,ur,sym)
            pxz(1,1)=0.0
            call vcfftb(pxz,pxz(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz,pxz(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call getxzp(pxz2,y,2,ur,sym)
            call vcfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            do z=1,nz
               do x=1,nx
                  pxz(x,z)=-pxz(x,z)*pxz2(x,z)/real(nx*nz)
               end do
            end do
            call vrfftf(pxz,pxz(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call vcfftf(pxz,pxz(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            if (ivar.eq.14) then
               call getxzp(pxz2,y,4,ur,sym)
               do z=1,nz
                  do x=1,nx
                     pxz(x,z)=pxz(x,z)+pxz2(x,z)/re
                  end do
               end do
            end if
            if (ivar.eq.15) then
               do z=1,nz
                  do x=1,nx
                     pxz(x,z)=pxz(x,z)*pxy(1,y)
                  end do
               end do
            end if
            call putxzp(pxz,y,4,ur)
         end do
      end if
c
c     Enstrophy <omega,omega>, modulus of vorticity = sqrt(enstrophy)
c
      if (ivar.eq.16.or.ivar.eq.17) then
c
c     en=wy*wy+vz*vz-2*wy*vz+uz*uz+wx*wx-2*uz*wx+vx*vx+uy*uy-2*vx*uy
c     step 1 calculate uy and wy in xy-plane put in 3-D ur(,,,4) and ur(,,,5)
c
         do i=1,2
            do z=1,nzc
               call getxyp(pxy,z,2*i-1,ur)
               call vchbf(pxy,w,nyp,nx,nx,1,prey)
               call rdcheb(pxy,nyp,nx,nx)
               call vchbb(pxy,w,nyp,nx,nx,1,prey)
               if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
                  do y=1,nyp
                     do x=1,nx/2
                        pxy(x,y)=pxy(x,y)*dstar*(2./real(nyp-1))
                     end do
                  end do
               else
                  do y=1,nyp
                     do x=1,nx/2
                        pxy(x,y)=pxy(x,y)*(2./real(nyp-1))
                     end do
                  end do
               end if
               call putxyp(pxy,z,3+i,ur)
            end do
         end do
c
c     Step 2 calculate enstrophy xz-plane by xz-plane
c     en=wy*wy+vz*vz-2*wy*vz+uz*uz+wx*wx-2*uz*wx+vx*vx+uy*uy-2*vx*uy
c
         do y=1,nyp
c
c     Get v,uy,wy
c
            sym=1.
            call getxzp(pxz,y,2,ur,sym)
            call getxzp(pxz2,y,4,ur,sym)
            sym=-1.
            call getxzp(pxz3,y,5,ur,sym)
c
c     Calculate omx=wy-vz and omz=vx-uy
c
            do z=1,nz
               do x=1,nx/2
                  pxz2(2*x-1,z)=pxz2(2*x-1,z)+beta(z)*pxz(2*x,z)
                  pxz2(2*x,z)=pxz2(2*x,z)-beta(z)*pxz(2*x-1,z)
                  pxz3(2*x-1,z)=-alfa(x)*pxz(2*x,z)-pxz3(2*x-1,z)
                  pxz3(2*x,z)=alfa(x)*pxz(2*x-1,z)-pxz3(2*x,z)
               end do
            end do
c
c     Transform into physical space
c
            call vcfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call vcfftb(pxz3,pxz3(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz3,pxz3(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
c
c     Accumulatie en=omy**2
c
            do z=1,nz
               do x=1,nx
                  pxz(x,z)=pxz2(x,z)**2+pxz3(x,z)**2
               end do
            end do
c
c     Get u,w
c
            sym=1.
            call getxzp(pxz2,y,1,ur,sym)
            sym=-1.
            call getxzp(pxz3,y,3,ur,sym)
c
c     Calculate omy=uz-wx
c
            do z=1,nz
               do x=1,nx/2
                  tr=alfa(x)*pxz3(2*x,z)-beta(z)*pxz2(2*x,z)
                  pxz2(2*x,z)=beta(z)*pxz2(2*x-1,z)-alfa(x)*
     &                 pxz3(2*x-1,z)
                  pxz2(2*x-1,z)=tr
               end do
            end do
c
c     Transform into physical space
c
            call vcfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call vrfftb(pxz2,pxz2(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            if (ivar.eq.16) then
c
c     Enstrophy
c
               do z=1,nz
                  do x=1,nx
                     pxz(x,z)=(pxz(x,z)+pxz2(x,z)**2)/real(nx*nz)
                  end do
               end do
            else
c
c     Modulus of vorticity
c
               do z=1,nz
                  do x=1,nx
                     pxz(x,z)=sqrt(pxz(x,z)+pxz2(x,z)**2)/real(nx*nz)
                  end do
               end do
            end if
            call vrfftf(pxz,pxz(2,1),wxz,wxz(2,1),nx,nz,2,nx+2,prex)
            call vcfftf(pxz,pxz(2,1),wxz,wxz(2,1),nz,nx/2,nx+2,2,prez)
            call putxzp(pxz,y,4,ur)
         end do
      end if
c
c     uw=u*cos(phiuw)+w*sin(phiuw)
c     If uwsub is true subtract laminar velocity
c
      if (ivar.eq.18) then
         do z=1,nzc
            call getxyp(pxy,z,1,ur)
            call getxyp(pxy2,z,3,ur)
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=cos(phiuw)*pxy(x,y)+sin(phiuw)*pxy2(x,y)
               end do
            end do
            if (z.eq.1.and.uwsub) then
               do y=1,nyp
                  if (fltype.eq.1.or.fltype.eq.4) then
                     pxy(1,y)=pxy(1,y)-(1.-eta(y)**2)-pxy(1,nyp)
                  else
                     pxy(1,y)=pxy(1,y)-eta(y)
                  end if
               end do
            end if
            call putxyp(pxy,z,4,ur)
         end do
      end if
c
c     Continuity dui/dxi
c
      if (ivar.eq.19) then
         do z=1,nzc
            call getxyp(pxy3,z,1,ur) ! get u
            call getxyp(pxy,z,2,ur)  ! get v
            call vchbf(pxy,w,nyp,nx,nx,1,prey)
            call rdcheb(pxy,nyp,nx,nx)
            if (fltype.lt.1.or.fltype.ge.6.or.fltype.eq.3) then
               do y=1,nyp
                  do x=1,nx/2
                     pxy(x,y)=pxy(x,y)*dstar
                  end do
               end do
            end if
            call vchbb(pxy,w,nyp,nx,nx,1,prey)
            call getxyp(pxy2,z,3,ur) ! get w
            do y=1,nyp
               do x=1,nx/2
                  pxy(x,y)=pxy(x,y)*(2./real(nyp-1))+
     &           (0.,1.)*alfa(x)*pxy3(x,y)+(0.,1.)*beta(z)*pxy2(x,y)
c            pxy(x,y)=(0.,1.)*alfa(x)*pxy3(x,y) ! ux
c            pxy(x,y)=pxy(x,y)*(2./real(nyp-1)) ! vy
c            pxy(x,y)=(0.,1.)*beta(z)*pxy2(x,y) ! wz
               end do
            end do
            call putxyp(pxy,z,4,ur)
         end do
      end if

      return

      end subroutine cvar
