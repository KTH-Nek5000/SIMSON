c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine stnois(ur,seed,nxn,nyn,nzn,ed,
     &     alfa,beta,eta,uxy)
c
c     Adds noise to a field stored in ur
c     the noise is composed of Stokes' modes with a total energy density
c     of ed evenly distributed among nxn*nyn*nzn*2 modes
c
c     There appears to be a problem with the divergence-free condition
c     of the noise based on Stokes modes. For high numbers of eigen-
c     functions (e.g. N>nyp/2) large errors in the divergence have
c     been observed. For lower numbers, however, no problems appear.
c
c     This issue is related to the computation of the
c     Stokes modes, which uses a potentially ill-conditioned algorithm.
c     Since there is no problems for low numbers of nyn, we do not
c     see the necessity to exchange/enhance the numerical computation
c     of the Stokes modes.
c

      implicit none

      include 'par.f'
c
c     Main storage
c
      complex ur(memnx,memny,memnz,3)
c
c     xy storage
c
      complex uxy(nx/2,nyp,3,2)
c
c     Parameters
c
      real ed
      integer nxn,nyn,nzn,seed
      real alfa(nx/2),beta(nz),eta(nyp)
c
c     Function
c
      real ran2
c
c     Local storage
c
      integer x,y,z,i,is,iz,zn,ni
      real wm,rnd,k2
      complex ph
      complex v(nyp,2),dv(nyp,2),omy(nyp,2)
      real vst(nyp),dvst(nyp),omyst(nyp),ust(nyp)
c
c     Local parameters
c
      complex im
      parameter (im=(0.,1.))
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Integration weights
c
      real wint(nyp)
c
c     Check input noise parameters
c
      if (nxn.lt.0.or.nxn.gt.nx/2) then
         write(*,*) 'Reduce nxn to maximum nx/2=',nx/2
         stop
      end if
      if (nzn.lt.0.or.nzn.gt.nz) then
         write(*,*) 'Reduce nzn to maximum nz=',nz
         stop
      end if
      if (nyn.lt.0.or.nyn.gt.nyp/2) then
         write(*,*) 'Warning:'
         write(*,*) 'nyn possibly too large to obtain'
         write(*,*) 'divergence-free noise.'
         write(*,*)
      end if


      write(*,*) 'Generating noise'
c
c     Set up integration weights
c
      do y=1,nyp
         wint(y)=-.5
         do ni=1,(nyp+1)/2-2
            wint(y) = wint(y) + cos(pi*real(2*ni*(y-1))/real(nyp-1))/
     &           real((2*ni)**2-1)
         end do
         wint(y) = wint(y) + .5*cos(pi*real(y-1))/real((nyp-1)**2-1)
         wint(y) = -4./real(nyp-1)*wint(y)
         if (y.eq.1.or.y.eq.nyp) wint(y)=wint(y)*.5
      end do
c
c     Calculate weight for each mode
c
      wm=sqrt(2.*ed/real((2*nxn-1)*nyn*nzn*2))
c
c     Start with beta= zero special case
c
      do i=1,3
         call getxyp(uxy(1,1,i,1),1,i,ur)
      end do
c
c     Special case wavenumber zero, stokes creates u/w eigenmodes
c
      do is=2,3
         do i=1,nyn/2
            call stokes(vst,dvst,ust,eta,wint,nyp,
     &           0.,0.,is,i)
            if (nfzsym.eq.1) then
c
c     For symmetry, w(0,0) must be zero so we are forced to randomize
c     the amplitude of the u-component. Thus the input energy density is not
c     satisified exactly, but rather is the expectancy value. Note the factor 2.,
c     one sqrt(2.) factor to compensate for that there is only one mode and one
c     because the the mean energy of the random direction projection is 1/2
c
               rnd=ran2(seed)
               do y=1,nyp
                  uxy(1,y,1,1)=uxy(1,y,1,1)+ust(y)*2.*wm*cos(2.*pi*rnd)
               end do
            else
c
c     No symmetry, let u and w be a constant length vector with random direction
c     the sqrt(2.) as above
c
               rnd=ran2(seed)
               do y=1,nyp
                  uxy(1,y,1,1)=uxy(1,y,1,1)+ust(y)*
     &                 sqrt(2.)*wm*cos(2.*pi*rnd)
                  uxy(1,y,3,1)=uxy(1,y,3,1)+ust(y)*
     &                 sqrt(2.)*wm*sin(2.*pi*rnd)
               end do
            end if
         end do
      end do
c
c     beta=0,alfa<>0
c
      do x=2,nxn
         do y=1,nyp
            v(y,1)=0.0
            dv(y,1)=0.0
            omy(y,1)=0.0
         end do
         do is=0,3
            do i=1,nyn/2
               call stokes(vst,dvst,omyst,eta,wint,nyp,
     &              alfa(x),0.,is,i)
               ph=exp(im*2.*pi*ran2(seed))
               if (is.le.1) then
c
c     Add v-mode with random phase
c
                  do y=1,nyp
                     v(y,1)=v(y,1)+vst(y)*wm*ph
                     dv(y,1)=dv(y,1)+dvst(y)*wm*ph
                  end do
               else
c
c     Add omy-mode with random phase
c
                  do y=1,nyp
                     omy(y,1)=omy(y,1)+omyst(y)*wm*ph
                  end do
               end if
            end do
         end do
c
c     Add this Fourier mode to the total field
c
         do y=1,nyp
c
c     For symmetry no energy is allowed in w for beta=0, then energy in u,v is
c     doubled
c
            if (nfzsym.eq.0) then
               uxy(x,y,1,1)=uxy(x,y,1,1)+im/alfa(x)*dv(y,1)
               uxy(x,y,2,1)=uxy(x,y,2,1)+v(y,1)
               uxy(x,y,3,1)=uxy(x,y,3,1)+im/alfa(x)*omy(y,1)
            else
               uxy(x,y,1,1)=uxy(x,y,1,1)+sqrt(2.)*im/alfa(x)*dv(y,1)
               uxy(x,y,2,1)=uxy(x,y,2,1)+sqrt(2.)*v(y,1)
            end if
         end do
      end do
      do i=1,3
         call putxyp(uxy(1,1,i,1),1,i,ur)
      end do
c
c     beta<>0 note treatment in pairs to ensure correct symmetry for alfa=0
c
      do z=2,nzn/2+1
         write(*,*) 'Making noise z= ', z
         zn=nz+2-z
         do i=1,3
            call getxyp(uxy(1,1,i,1),z,i,ur)
            if (nfzsym.eq.0) call getxyp(uxy(1,1,i,2),zn,i,ur)
         end do
         do x=1,nxn
            k2=alfa(x)**2+beta(z)**2
            do iz=1,2-nfzsym
               do y=1,nyp
                  v(y,iz)=0.0
                  dv(y,iz)=0.0
                  omy(y,iz)=0.0
               end do
            end do
c
c     Generate the Stokes modes
c
            do is=0,3
               do i=1,nyn/2
                  call stokes(vst,dvst,omyst,eta,wint,nyp,
     &                 alfa(x),beta(z),is,i)
                  do iz=1,2-nfzsym
                     ph=exp(im*2.*pi*ran2(seed))
                     if (is.le.1) then
c
c     Add v-mode with random phase
c
                        do y=1,nyp
                           v(y,iz)=v(y,iz)+vst(y)*wm*ph
                           dv(y,iz)=dv(y,iz)+dvst(y)*wm*ph
                        end do
                     else
c
c     Add omy-mode with random phase
c
                        do y=1,nyp
                           omy(y,iz)=omy(y,iz)+omyst(y)*wm*ph
                        end do
                     end if
                  end do
               end do
            end do
c
c     Add this Fourier mode to the total field
c
            if (x.eq.1) then
c
c     Special case alfa=0, note that the wavenumbers with opposite beta
c     must have complex conjugated coefficients for the velocity to be real
c     for the symmetric case the waves must also hold their symmetry
c     this means that u,v should be real and w imaginary
c     this again cannot be fulfilled with exact energy, so that
c     energy is only correct for an expectancy value
c
               if (nfzsym.eq.0) then
                  do y=1,nyp
                     uxy(x,y,1,1)=uxy(x,y,1,1)-im/k2*beta(z)*omy(y,1)
                     uxy(x,y,2,1)=uxy(x,y,2,1)+v(y,1)
                     uxy(x,y,3,1)=uxy(x,y,3,1)+im/k2*beta(z)*dv(y,1)
                     uxy(x,y,1,2)=uxy(x,y,1,2)-conjg(im/k2*beta(z)*
     &                            omy(y,1))
                     uxy(x,y,2,2)=uxy(x,y,2,2)+conjg(v(y,1))
                     uxy(x,y,3,2)=uxy(x,y,3,2)+conjg(im/k2*beta(z)*
     &                            dv(y,1))
                  end do
               else
                  do y=1,nyp
                     uxy(x,y,1,1)=uxy(x,y,1,1)-
     &                    real(im/k2*beta(z)*omy(y,1))*sqrt(2.)
                     uxy(x,y,2,1)=uxy(x,y,2,1)+real(v(y,1))*sqrt(2.)
                     uxy(x,y,3,1)=uxy(x,y,3,1)+
     &                    im*real(-im*(im/k2*beta(z)*dv(y,1)))*sqrt(2.)
                  end do
               end if
            else
c
c     alfa<>0 no coupling for positive and negative beta
c
               do y=1,nyp
                  uxy(x,y,1,1)=uxy(x,y,1,1)+
     &                 im/k2*(alfa(x)*dv(y,1)-beta(z)*omy(y,1))
                  uxy(x,y,2,1)=uxy(x,y,2,1)+v(y,1)
                  uxy(x,y,3,1)=uxy(x,y,3,1)+
     &                 im/k2*(beta(z)*dv(y,1)+alfa(x)*omy(y,1))
               end do
               if (nfzsym.eq.0) then
                  do y=1,nyp
                     uxy(x,y,1,2)=uxy(x,y,1,2)+
     &                    im/k2*(alfa(x)*dv(y,2)-beta(zn)*omy(y,2))
                     uxy(x,y,2,2)=uxy(x,y,2,2)+v(y,2)
                     uxy(x,y,3,2)=uxy(x,y,3,2)+
     &                    im/k2*(beta(zn)*dv(y,2)+alfa(x)*omy(y,2))
                  end do
               end if
            end if
         end do
         do i=1,3
            call putxyp(uxy(1,1,i,1),z,i,ur)
            if (nfzsym.eq.0) call putxyp(uxy(1,1,i,2),zn,i,ur)
         end do
      end do
      write(*,*)

      end subroutine stnois
