c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine steadywallrough(bu1,bu2,wallur,wallui,
     &     wallwr,wallwi,v_wall,wlvr,wlvi,prexn,prey,prezn,presn,prean,
     &     alfa,beta,zl,xsc,zsc,dstar,xst,xlen,delx,h_rough,xstart,xend,
     &     xrise,xfall,zrand,hrms,nzr,rseed,rzd,psi,tpsi,taylor4,rghfil,
     &     roughfile,my_node)
c
c     Calculates the boundary condition at y=0 resulting from
c     surface roughness. The latter is modelled in terms of a Taylor
c     series (projection of the no-slip conditions from the desired
c     bump contour to the "wall" y=0)
c
c     Taylor series is truncated at 1st order or at 4th order (if
c     taylor4==T)
c
c     The projected wall-boundary conditions are calculated using
c     the BASE FLOW (bu1,bu2). Hence, this routine is called only once
c     before the time-step loop begins.
c
c     CALLED by bla.f (if wbci==-1)
c
c     CALLS vchbf, vchbb,dcheb,rdcheb,fft2df
c
      implicit none

      include 'par.f'

#ifdef MPI
      include 'mpif.h'
#endif

c
c     Parameters
c
      integer nxm,nxn,nzn,lenx
      parameter(nxm=nxp/2+1,nxn=nx/2+1,nzn=nz/2+1)
      parameter (lenx=nx/8)

      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Global variables
c
      real bu1(nxm,nyp,3),bu2(nxm,nyp,3)

      integer xst,xlen,nzr,rseed,my_node

      real zl,xsc,zsc,psi,tpsi,dstar,delx
      real h_rough,xstart,xend,xrise,xfall

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15)
      real prean(nz*3/4+15),prey(nyp*2+15),w(nxm,nyp)
      real workr(nz+15)

      real alfa(nx/2*mbz),beta(nz)

      real rzd(nz+2)

      logical v_wall,taylor4,rghfil,hrms,zrand
      character(32) roughfile
c
c     Output
c
      real wallur(nxm,memnz),wallui(nxm,memnz)
      real wlvr(nxm,memnz),wlvi(nxm,memnz)
      real wallwr(nxm,memnz),wallwi(nxm,memnz)
c
c     Local variables
c
      integer x,xsh,xind,y,z,zb,i,split,len

      real wr(nxn*mby,nz),wi(nxn*mby,nz),wz(nz+2)

      real lwallur(nxn*mby,nz),lwallui(nxn*mby,nz)
      real lwallvr(nxn*mby,nz),lwallvi(nxn*mby,nz)
      real lwallwr(nxn*mby,nz),lwallwi(nxn*mby,nz)

      real hfunc(lenx,nz,2)
      real arg,hr
      real xc11,xc22,xco,fx1,fx2,umx,vmx,wmx
      real step,tmpr1,tmpr2
c      real zshift1,zshift2
      real delty,norm,hd

      real app1(nyp),app2(nyp),app21(nyp),app22(nyp)
      real app31(nyp),app32(nyp),app41(nyp),app42(nyp)
      real du1wall(lenx,3),du2wall(lenx,3)
      real d2u1wall(lenx,3),d2u2wall(lenx,3)
      real d3u1wall(lenx,3),d3u2wall(lenx,3)
      real d4u1wall(lenx,3),d4u2wall(lenx,3)

      complex ui

      character(36) rghfile


      if (my_node.eq.0) then
         write(ios,*)'--> Wall-rough. BCs set (based on initial flow).'
         if (psi.eq.90) then
            write(ioe,*)'! Please set psi < 90 degrees !'
            stop
         end if
      end if
c
c Imaginary unit
c
      ui=(0.,1.)
c
c     Angle by which roughness crests are turned
c
      psi  = pi*psi/180
      tpsi = tan(psi)
c
c     Initialize the wall B.C.s:
c
      lwallur(:,:)=0.
      lwallui(:,:)=0.
      lwallvr(:,:)=0.
      lwallvi(:,:)=0.
      lwallwr(:,:)=0.
      lwallwi(:,:)=0.
c
c     Calculation of the vel. derivatives at y=0
c
      delty=2./real(nyp-1)      ! y step width (internal scaling)

      xsh = xst+nxp/4

      do i=1,3
         do x=1,xlen
            xind=x+xsh-1

            do y=1,nyp
               app1(y)=bu1(xind,y,i)
               app2(y)=bu2(xind,y,i)
            end do

            call vchbf(app1,w,nyp,1,1,1,prey)
            call vchbf(app2,w,nyp,1,1,1,prey)
c
c     Normalise
c
            do y=1,nyp
               app1(y)=app1(y)*delty
               app2(y)=app2(y)*delty
            end do

            call rdcheb(app1,nyp,1,1) !1st deriv. of vel. profs.
            call rdcheb(app2,nyp,1,1)

            if (taylor4) then

               call dcheb(app21,app1,nyp,1,1) ! 2nd derivative
               call dcheb(app22,app2,nyp,1,1)

               call dcheb(app31,app21,nyp,1,1) ! 3rd derivative
               call dcheb(app32,app22,nyp,1,1)

               call dcheb(app41,app31,nyp,1,1) ! 4th derivative
               call dcheb(app42,app32,nyp,1,1)

               call vchbb(app21,w,nyp,1,1,1,prey)
               call vchbb(app22,w,nyp,1,1,1,prey)

               call vchbb(app31,w,nyp,1,1,1,prey)
               call vchbb(app32,w,nyp,1,1,1,prey)

               call vchbb(app41,w,nyp,1,1,1,prey)
               call vchbb(app42,w,nyp,1,1,1,prey)

               d2u1wall(x,i)=app21(nyp) !2nd derivative
               d2u2wall(x,i)=app22(nyp)

               d3u1wall(x,i)=app31(nyp) !3rd derivative
               d3u2wall(x,i)=app32(nyp)

               d4u1wall(x,i)=app41(nyp) !4th derivative
               d4u2wall(x,i)=app42(nyp)
            end if

            call vchbb(app1,w,nyp,1,1,1,prey)
            call vchbb(app2,w,nyp,1,1,1,prey)

            du1wall(x,i)=app1(nyp) !1st derivative
            du2wall(x,i)=app2(nyp)
         end do
      end do
c
c     Bump function and projected BCs
c
      xsh = xst-0.5*nxn
      if (xsh.le.0) xsh = xsh+nxn
c
c     Check whether proper # of modes
c
      if (my_node.eq.0) then
         if (nzr.lt.0) then
            write(ioe,*) 'Random roughn.: nzr must be .ge. 0: ',nzr
            stop
         end if
         if (nzr.gt.nz/2) then
            write(ioe,*) 'Random roughn.: nzr must be .le. nz/2: ',
     &           nz/2,nzr
            stop
         end if
      end if
c
c     Generate a random z distribution
c
      rzd(:) = 0.
      nzr = 2*nzr
      if (nzr.eq.0) then
         nzr = 1
      else
         call zdran(rzd,hrms,zrand,nzr,rseed,workr)
c
c     Fourier shift the spanw. dependence to zsc
c
         if (nz.gt.1.and.nfzsym.eq.0) then
            if (zsc.ne.0.) then
               call vrfftf(rzd,rzd(2),wz,wz(2),nz,1,2,1,workr)
               do z=1,nzr/2+1
                  arg=zsc/zl*2.*pi*real(z-1)
                  hr=rzd(2*z-1)*cos(arg)-rzd(2*z)*sin(arg)
                  rzd(2*z)=rzd(2*z)*cos(arg)+rzd(2*z-1)*sin(arg)
                  rzd(2*z-1)=hr
               end do
               call vrfftb(rzd,rzd(2),wz,wz(2),nz,1,2,1,workr)
               do z=1,nz
                  rzd(z)=1./real(nz)*rzd(z)
               end do
            end if
         end if

      end if


      do z=1,nz
         if (nzr.eq.1) then
            tmpr1 = 1.
            tmpr2 = 1.
         else
            tmpr1 = rzd(z)
            tmpr2 = tmpr1
         end if

         do x=1,xlen
            xind=x+xsh-1

            xc11=xstart+(2*x-2)*delx+xsc
            fx1=step((xc11-xstart)/xrise)-step((xc11-xend)/xfall+1)

            xc22 = xstart+(2*x-1)*delx+xsc
            fx2=step((xc22-xstart)/xrise)-step((xc22-xend)/xfall+1)

            hfunc(x,z,1)=h_rough*fx1*tmpr1
            hfunc(x,z,2)=h_rough*fx2*tmpr2

            lwallur(xind,z)=-hfunc(x,z,1)*du1wall(x,1) !"roughness
            lwallui(xind,z)=-hfunc(x,z,2)*du2wall(x,1) !projection"

            lwallwr(xind,z)=-hfunc(x,z,1)*du1wall(x,3)
            lwallwi(xind,z)=-hfunc(x,z,2)*du2wall(x,3)

            if (taylor4) then !Take 2nd, 3rd, 4th deriv. into account
               lwallur(xind,z)=lwallur(xind,z)-
     &              0.5*(hfunc(x,z,1)**2)*d2u1wall(x,1)-
     &              (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,1)-
     &              (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,1)
               lwallui(xind,z)=lwallui(xind,z)-
     &              0.5*(hfunc(x,z,2)**2)*d2u2wall(x,1)-
     &              (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,1)-
     &              (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,1)

               lwallwr(xind,z)=lwallwr(xind,z)-
     &              0.5*(hfunc(x,z,1)**2)*d2u1wall(x,3)-
     &              (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,3)-
     &              (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,3)
               lwallwi(xind,z)=lwallwi(xind,z)-
     &              0.5*(hfunc(x,z,2)**2)*d2u2wall(x,3)-
     &              (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,3)-
     &              (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,3)
            end if

            if (v_wall) then
               lwallvr(xind,z)=-hfunc(x,z,1)*du1wall(x,2)
               lwallvi(xind,z)=-hfunc(x,z,2)*du2wall(x,2)
               if (taylor4) then !Take 2nd, 3rd, 4th deriv. into account
                  lwallvr(xind,z)=lwallvr(xind,z)-
     &                 0.5*(hfunc(x,z,1)**2)*d2u1wall(x,2)-
     &                 (1./6.)*(hfunc(x,z,1)**3)*d3u1wall(x,2)-
     &                 (1./24.)*(hfunc(x,z,1)**4)*d4u1wall(x,2)
                  lwallvi(xind,z)=lwallvi(xind,z)-
     &                 0.5*(hfunc(x,z,2)**2)*d2u2wall(x,2)-
     &                 (1./6.)*(hfunc(x,z,2)**3)*d3u2wall(x,2)-
     &                 (1./24.)*(hfunc(x,z,2)**4)*d4u2wall(x,2)
               end if
             end if
          end do
       end do
c
c     Write u,v,w values (at wall, where bump height max.) on screen:
c
       if (my_node .eq. 0) then
          x = 0.5*(2*xsh+xlen-0.25)

          umx = 0.
          vmx = 0.
          wmx = 0.

          do z=1,nz
             umx=min(umx,lwallur(x,z))
             vmx=min(vmx,lwallvr(x,z))
             wmx=min(wmx,lwallwr(x,z))
          end do

          write(ios,*)'--> u,v,w|0 at h_max:'
          write(ios,118) umx,vmx,wmx
c
c     Write the roughness function into a file?
c
          if (rghfil) then

             split = nz/4

             hd=0.

             call CharStringLen(roughfile,len)
             call NewCharString(roughfile,len,'xxx',rghfile)

             open(unit=121,status='replace',file=rghfile)

             call NewCharString(roughfile,len,'hhh',rghfile)

             open(unit=122,status='replace',file=rghfile,
     &            form='formatted')

             do x=1,xst-1
                xc11=(2*x-2)*delx+xsc
                xc22=(2*x-1)*delx+xsc
                write(121,*) xc11/dstar
                write(121,*) xc22/dstar
                do i = 1,split
                   write(122,119) (hd, z= (i-1)*4+1,i*4)
                   write(122,119) (hd, z= (i-1)*4+1,i*4)
                end do
             end do

             xco = xc22

             do x=1,xlen
                xc11=xco+(2*x-1)*delx+xsc
                xc22=xco+2*x*delx+xsc
                write(121,*) xc11/dstar
                write(121,*) xc22/dstar
                do i = 1,split
                   write(122,119) (hfunc(x,z,1)/dstar, z= (i-1)*4+1,i*4)
                end do
                do i = 1,split
                   write(122,119) (hfunc(x,z,2)/dstar, z= (i-1)*4+1,i*4)
                end do
             end do

             do x=xst+xlen,nxn-1
                xc11=(2*x-2)*delx+xsc
                xc22=(2*x-1)*delx+xsc
                write(121,*) xc11/dstar
                write(121,*) xc22/dstar
                do i = 1,split
                   write(122,119) (hd, z= (i-1)*4+1,i*4)
                   write(122,119) (hd, z= (i-1)*4+1,i*4)
                end do
             end do

             close(121)
             close(122)
          end if
       end if

 118   format(3F22.16)
 119   format(4F16.8)
c
c     Transform to Fourier space
c
      call sfft2df(lwallur,lwallui,.true.,1,prexn,prezn,presn,prean,
     &      wr,wi)
      call sfft2df(lwallwr,lwallwi,.false.,1,prexn,prezn,presn,prean,
     &     wr,wi)
      if (v_wall) call sfft2df(lwallvr,lwallvi,.true.,1,
     &     prexn,prezn,presn,prean,wr,wi)
c
c     Normalization
c
      norm=real(nx*nz)
      do z=1,nz/2
         do x=1,nx/2
            lwallur(x,z)=lwallur(x,z)/norm
            lwallui(x,z)=lwallui(x,z)/norm
            lwallwr(x,z)=lwallwr(x,z)/norm
            lwallwi(x,z)=lwallwi(x,z)/norm
         end do
      end do

      if (nfzsym.eq.0) then
         do z=nz/2+1,nz
            do x=1,nx/2
               lwallur(x,z)=lwallur(x,z)/norm
               lwallui(x,z)=lwallui(x,z)/norm
               lwallwr(x,z)=lwallwr(x,z)/norm
               lwallwi(x,z)=lwallwi(x,z)/norm
            end do
         end do
      end if
c
c     If bump condition for v is also prescribed
c
      if (v_wall) then
         do z=1,nz/2
            do x=1,nx/2
               lwallvr(x,z)=lwallvr(x,z)/norm
               lwallvi(x,z)=lwallvi(x,z)/norm
            end do
         end do

         if (nfzsym.eq.0) then
            do z=nz/2+1,nz
               do x=1,nx/2
                  lwallvr(x,z)=lwallvr(x,z)/norm
                  lwallvi(x,z)=lwallvi(x,z)/norm
               end do
            end do
         end if
      end if
c
c     From this point the subroutine is parallelized
c     Obtain Dv and eta
c
      do z=1,memnz
         zb=my_node*memnz+z
c
c     Re{Dv},Im{Dv},Re{eta} and Im{eta}
c
         do x=1,nx/2
            wallur(x,z)=alfa(x)*lwallui(x,zb)+beta(zb)*lwallwi(x,zb)
            wallui(x,z)=-alfa(x)*lwallur(x,zb)-beta(zb)*lwallwr(x,zb)
            wallwr(x,z)=alfa(x)*lwallwi(x,zb)-beta(zb)*lwallui(x,zb)
            wallwi(x,z)=-alfa(x)*lwallwr(x,zb)+beta(zb)*lwallur(x,zb)
         end do
      end do
c
c     Re{v},Im{v}
c
      if (v_wall) then
         do z=1,memnz
            zb=my_node*memnz+z
            do x=1,nx/2
               wlvr(x,z)=lwallvr(x,zb)
               wlvi(x,z)=lwallvi(x,zb)
            end do
         end do
      end if
c
c     In (1,1), i.e. (alfa,beta)=(0,0) the values of u and w
c     and not Dv and eta are put back (not necessary for v)
c
      if (my_node.eq.0) then
         wallur(1,1)=lwallur(1,1)
         wallwr(1,1)=lwallwr(1,1)
      end if


      end subroutine steadywallrough


!_________________________________________________________


      subroutine zdran(ttz,hrms,zrand,nzt,seed,workr)
c
c     Create a function ttz with random dependence of z
c     containing nzt Fourier components
c
      implicit none

      include 'par.f'

      integer nzt,seed
      real ttz(nz+2),workr(nz+15)
      logical hrms,zrand
c
c     Local variables
c
      integer z
      real norm,tmean,trms,ttzmx
      complex ph
      real w(nz+2)
c
c     Local parameters
c
      integer zend
      complex im
      parameter (im=(0.,1.))
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     External function
c
      real ran2
c
c     Initialise function and work array
c
      ttz(:)=0.0
      call vrffti(nz,workr,0)
      zend = nzt/2+1
c
c     Construct ttz with alternating real/imaginary parts
c     up to nzt/2+1 and zeros at positions 2 and nz+2
c
      if (nfzsym.eq.0) then
c
c     For the non-symmetric case, we take a constant
c     amplitude and randomize the phase
c
         if (zrand) then
            do z=1,zend
               ph = exp(im*2.*pi*ran2(seed))
               ttz(2*z-1) = real(ph)
               ttz(2*z) = aimag(ph)
            end do
         else
            ph = exp(im*0.5*pi)
            ttz(2*zend-1) = real(ph)
            ttz(2*zend) = aimag(ph)
          end if
          ttz(nz+2) = 0.0
      else
c
c     For the symmetric case the transform should be real
c     which leads to random amplitude instead of random phase
c     we raise the amplitude to compensate
c
         if (zrand) then
            do z=1,zend
               ph = exp(im*2.*pi*ran2(seed))
               ttz(2*z-1) = real(ph)*sqrt(2.)
               ttz(2*z) = 0.0
            end do
         else
            ph = exp(im*0.5*pi)
            ttz(2*zend-1) = real(ph)*sqrt(2.)
            ttz(2*zend) = 0.0
         end if
      end if
c
c     Transform to real space
c
      if (nz.gt.1) call vrfftb(ttz(1),ttz(2),w(1),w(2),nz,1,2,1,workr)
c
c     Compute/subtract mean and compute max of ttz:
c
      tmean = 0.
      do z = 1,nz
         tmean = tmean+ttz(z)
      end do
      tmean = tmean/real(nz)
      if (nzt.gt.1) then
         do z = 1,nz
            ttz(z) = ttz(z)-tmean
         end do
      end if
      ttzmx = 0.
      do z = 1,nz
         ttzmx = max(abs(ttz(z)),ttzmx)
      end do
c
c     Normalize and compute rms of ttz:
c
      norm = 1.0/ttzmx
      do z = 1,nz
         ttz(z) = norm*ttz(z)
      end do
      if (hrms) then
         trms = 0.
         do z = 1,nz
            trms = trms+ttz(z)**2
         end do
         trms = trms/real(nz)
         trms = sqrt(trms)
         norm = 1.0/trms
         do z = 1,nz
            ttz(z) = norm*ttz(z)
         end do
      end if


      end subroutine zdran


!_____________________________________________________


! 2 routines creating new filenames from a given name

      subroutine CharStringLen(namin,len)

      implicit none

      character(32),intent(in) :: namin
      integer(4),intent(out) :: len

      ! Local:
      integer(4) :: i

      i=32
      do while (namin(i:i) .eq. ' ')
         i=i-1
      end do
      len=i

      end subroutine  CharStringLen

!_____________________________________________________

      subroutine NewCharString(namin,len,ext,namout)

      implicit none

      integer(4),intent(in) :: len

      character(32),intent(in) :: namin
      character(3),  intent(in) :: ext

      character(36),intent(out) :: namout

      !local:
      character(len) :: naminnew

      naminnew=namin
      namout=naminnew//'.'//ext

      end subroutine NewCharString
