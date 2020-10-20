c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine init_filter(cutoff,nx,ny,nz,xkoord,
     &     gewy1,gewy2,gewspecx,gewspecz,diags,my_node)
c
c     Initialise the LES filters based on the implementation by
c     Steffen Stolz, see Stolz (PhD thesis, 2000) and
c     Stolz, Adams, Kleiser (Phys. Fluids, 2001).
c
      implicit none

      real,    intent(in)  :: cutoff
      integer, intent(in)  :: nx,ny,nz
      real,    intent(in)  :: xkoord(ny)
      real,    intent(out) :: gewy1(ny,5),gewy2(ny,5)
      real,    intent(out) :: gewspecx(nx+1),gewspecz(nz+1)
      real,    intent(out) :: diags(ny,5)

      real gew1hom(5), gew2hom(5)

      real dmmm,dmm,dm,dp,dpp,dppp
      real ammx,amx,anx,apx,appx
      real amm,am,an,ap,app

      real gew1tmp(5), gew2tmp(5)
      real pi
      real rreGimpl,rimGimpl,rk,ra,rb,rc,rd
      real fact3z,factm3z

      integer i,k,l,ifail,my_node,ioutput

      ioutput = 0

      pi = 4.*atan(1.)
      gewspecx=0.
      gewspecz=0.

      if (my_node.eq.0) then
         write(*,*)
         write(*,*)
         write(*,*) '>>>>>>  Initialize primary LES filter  <<<<<<<'
         write(*,*) '--------------------------------------------'//
     &              '-----------------------'
         write(*,*) '* nx = ',nx
         write(*,*) '* ny = ',ny
         write(*,*) '* nz = ',nz
      end if

      if (mod(ny,2).eq.0) then
c
c     Although it would work
c
         write(*,*) 'ny is even. Stop.'
         stop
      end if
c
c     Primary filter for the homogeneous directions
c
c      write(*,*) '* cutoff = pi * ',cutoff
      if (cutoff.eq.0.) then
         write(*,*) 'Stop.'
      end if

      dmm = 2.
      dm  = 1.
      dp  = 1.
      dpp = 2.
      call coeff_filter_centered(dmm,dm,dp,dpp,cutoff,gew1hom,gew2hom,0)

      ! explicit filter
      ammx = gew1hom(1)
      amx  = gew1hom(2)
      anx  = gew1hom(3)
      apx  = gew1hom(4)
      appx = gew1hom(5)

      ! implicit filter
      amm  = gew2hom(1)
      am   = gew2hom(2)
      an   = gew2hom(3)
      ap   = gew2hom(4)
      app  = gew2hom(5)

      ! some writing
      if (my_node.eq.0) then
         write(*,'(a4,a6,9a10)') ' ','Y','M0','M1','M2','M3','M4',
     &        'SX'
         write(*,'(a4,a6,9e10.2)') 'h ',' ',
     &        (ammx*(-dmm)**l+amx*(-dm)**l+anx*(0.)**l+
     &        apx*(dp)**l+appx*(dpp)**l,l=0,4),
     &        ammx+amx+anx+apx+appx
         if (ioutput.ne.0) then
            call write_filter_centered( gew1hom, gew2hom,62 )
         end if
      end if

c---  convert filter into Fourier space
      do i=1,nx
         rk = (i-1)*pi/(nx-1)

         ! Fourier transform of the implicit filter
         ra = (ammx*cos(-rk*2.)+amx*cos(-rk)+
     &        anx+apx*cos(rk)+appx*cos(rk*2.))
         rc = (amm*cos(-rk*2.)+am*cos(-rk)+
     &        an+ap*cos(rk)+app*cos(rk*2.))
         rb = (ammx*sin(-rk*2.)+amx*sin(-rk)+
     &        apx*sin(rk)+appx*sin(rk*2.))
         rd = (amm*sin(-rk*2.)+am*sin(-rk)+
     &        ap*sin(rk)+app*sin(rk*2.))

         rreGimpl = (ra*rc+rb*rd)/(rc*rc+rd*rd)
         rimGimpl = (rb*rc-ra*rd)/(rc*rc+rd*rd)
         ! we know, however, that the imaginary part is zero
         rimGimpl = 0.

         gewspecx(i) = rreGimpl

      end do

      do i=1,nz
         rk = (i-1)*pi/(nz-1)

         ! Fourier transform of the implicit filter
         ra = (ammx*cos(-rk*2.)+amx*cos(-rk)+
     &        anx+apx*cos(rk)+appx*cos(rk*2.))
         rc = (amm*cos(-rk*2.)+am*cos(-rk)+
     &        an+ap*cos(rk)+app*cos(rk*2.))
         rb = (ammx*sin(-rk*2.)+amx*sin(-rk)+
     &        apx*sin(rk)+appx*sin(rk*2.))
         rd = (amm*sin(-rk*2.)+am*sin(-rk)+
     &        ap*sin(rk)+app*sin(rk*2.))

         rreGimpl = (ra*rc+rb*rd)/(rc*rc+rd*rd)
         rimGimpl = (rb*rc-ra*rd)/(rc*rc+rd*rd)
         ! we know, however, that the imaginary part is zero
         rimGimpl = 0.

         gewspecz(i) = rreGimpl

      end do


c---- and now for the inhomogeneous direction z


c     fill boundary points k=1
      k=1
      do l = 1, 5
         gewy1(1,l) = 0.
         gewy2(1,l) = 0.
      end do
      gewy1(1,1) = 1.
      gewy2(1,1) = 1.

      if (my_node.eq.0) then
         write(*,'(i4,f6.2,a)') k,xkoord(k),'   *** no filter ***'
      end if

      ! asymmetric filter from 1 to 5 with center 2
      ! for k=2
      k = 2
      ! get coefficients
      dm = xkoord(k-1) - xkoord(k)
      dp = xkoord(k)   - xkoord(k+1)
      dpp = xkoord(k)   - xkoord(k+2)
      dppp = xkoord(k)   - xkoord(k+3)
      call coeff_filter_r1_l3(dm,dp,dpp,dppp,cutoff,gew1tmp,gew2tmp,1)

      ! and put at position 2
      do l = 1, 5
         gewy1(k,l)      = gew1tmp(l)
         gewy2(k,l)      = gew2tmp(l)
      end do

      ! some writing
      if (my_node.eq.0) then
         write(*,'(i4,f6.2,9e10.2)') k,xkoord(k),
     &        (gew1tmp(1)*(-dm)**l+gew1tmp(2)*0.**l+gew1tmp(3)*(dp)**l+
     &        gew1tmp(4)*(dpp)**l  +gew1tmp(5)*(dppp)**l,l=0,4),
     &        gew1tmp(1)+gew1tmp(2)+gew1tmp(3)+gew1tmp(4)+gew1tmp(5)
         if (ioutput.ne.0) then
            call write_filter_r1_l3(gew1tmp,gew2tmp,64)
         end if
      end if



      ! and filter weights for the inner coordinates
      ! centered filter
      ! (for symmetry it could be improved)
c      do k = ny/2+1, ny-2
      do k = 3, ny-2
         ! get positions (always positive)
         dmm = xkoord(k-2) - xkoord(k)
         dm  = xkoord(k-1) - xkoord(k)
         dp  = xkoord(k)   - xkoord(k+1)
         dpp = xkoord(k)   - xkoord(k+2)

         call coeff_filter_centered(dmm,dm,dp,dpp,cutoff,gew1tmp,
     &        gew2tmp,0)

         ! and put at position k and N3N+2-k
         do l = 1, 5
            gewy1(k,l)          = gew1tmp(l)
            gewy2(k,l)          = gew2tmp(l)
c     for symmetry
c            gewy1(ny+1-k,6-l)   = gew1tmp(l)
c            gewy2(ny+1-k,6-l)   = gew2tmp(l)
         end do

         ! explicit filter
         ammx = gew1tmp(1)
         amx  = gew1tmp(2)
         anx  = gew1tmp(3)
         apx  = gew1tmp(4)
         appx = gew1tmp(5)

         ! implicit filter
         amm  = gew2tmp(1)
         am   = gew2tmp(2)
         an   = gew2tmp(3)
         ap   = gew2tmp(4)
         app  = gew2tmp(5)

         ! some writing
         if (my_node.eq.0) then
            write(*,'(i4,f6.2,9e10.2)') k,xkoord(k),
     &           (ammx*(-dmm)**l+amx*(-dm)**l+anx*(0.)**l+
     &           apx*(dp)**l+appx*(dpp)**l,l=0,4),
     &           ammx+amx+anx+apx+appx

            if (ioutput.ne.0) then
               call write_filter_centered(gew1tmp,gew2tmp,63)
            end if
         end if
      end do

      ! asymmetric filter from 1 to 5 with center 2
      ! for k=ny-1
      k = ny-1
      ! get coefficients
      dmmm = xkoord(k-3) - xkoord(k)
      dmm  = xkoord(k-2) - xkoord(k)
      dm   = xkoord(k-1) - xkoord(k)
      dp   = xkoord(k)   - xkoord(k+1)
      call coeff_filter_r3_l1(dmmm,dmm,dm,dp,cutoff,gew1tmp,gew2tmp,1)

      ! and put at position ny-1
      do l = 1, 5
         gewy1(k,l) = gew1tmp(l)
         gewy2(k,l) = gew2tmp(l)
      end do

      ! some writing
      if (my_node.eq.0) then
         write(*,'(i4,f6.2,9e10.2)') k,xkoord(k),
     &        (gew1tmp(1)*(-dmmm)**l+gew1tmp(2)*(-dmm)**l+
     &        gew1tmp(3)*(-dm)**l+gew1tmp(4)*(0)**l+
     &        gew1tmp(5)*(dp)**l,l=0,4),
     &        gew1tmp(1)+gew1tmp(2)+gew1tmp(3)+gew1tmp(4)+gew1tmp(5)
         if (ioutput.ne.0) then
            call write_filter_r3_l1(gew1tmp,gew2tmp,64)
         end if
      end if

c     fill boundary points k=ny
      k=ny
      do l = 1, 5
         gewy1(ny,l) = 0.
         gewy2(ny,l) = 0.
      end do
      gewy1(ny,5) = 1.
      gewy2(ny,5) = 1.
      if (my_node.eq.0) then
         write(*,'(i4,f6.2,a)') k,xkoord(k),'   *** no filter ***'
      end if


      if (my_node.eq.0.and.ioutput.ne.0) then
         ! output of filter weights
         ! homogeneous filter
         write(61,*) 'Filter weights:'
         write(61,*) 'Homogeneous filter'
         write(61,'(5f15.7)') (gew1hom(l),l=1,5)
         write(61,*) '&'
         write(61,'(5f15.7)') (gew2hom(l),l=1,5)
         write(61,*)
         write(61,*) 'Inhomogeneous filter'
         ! inhomogeneous filter
         do k = 1, ny
            write(61,'(i5,5F15.7)') k, (gewy1(k,l), l = 1, 5)
         end do
         write(61,*) '&'
         do k = 1, ny
            write(61,'(i5,5F15.7)') k, (gewy2(k,l) , l = 1, 5)
         end do
      end if
c
c     Compose diagonals and LU decomposition
c
      ! boundary points
      diags(1,1) = 0.
      diags(1,2) = 0.
      diags(1,3) = 1.
      diags(1,4) = 0.
      diags(1,5) = 0.

      diags(ny,1) = 0.
      diags(ny,2) = 0.
      diags(ny,3) = 1.
      diags(ny,4) = 0.
      diags(ny,5) = 0.

      ! asymmetric stencil (k=2 and k=ny-1)
      fact3z  = -gewy2(2,5)    / gewy2(3,5)
      factm3z = -gewy2(ny-1,1) / gewy2(ny-2,1)

      diags(2,1) = 0.
      diags(2,2) = gewy2(2,1) + gewy2(3,1)*fact3z
      diags(2,3) = gewy2(2,2) + gewy2(3,2)*fact3z
      diags(2,4) = gewy2(2,3) + gewy2(3,3)*fact3z
      diags(2,5) = gewy2(2,4) + gewy2(3,4)*fact3z

      diags(ny-1,1) = gewy2(ny-1,2) + gewy2(ny-2,2)*factm3z
      diags(ny-1,2) = gewy2(ny-1,3) + gewy2(ny-2,3)*factm3z
      diags(ny-1,3) = gewy2(ny-1,4) + gewy2(ny-2,4)*factm3z
      diags(ny-1,4) = gewy2(ny-1,5) + gewy2(ny-2,5)*factm3z
      diags(ny-1,5) = 0.

      ! symmetric stencil
      do k=3,ny-2
         diags(k,1) = gewy2(k,1)
         diags(k,2) = gewy2(k,2)
         diags(k,3) = gewy2(k,3)
         diags(k,4) = gewy2(k,4)
         diags(k,5) = gewy2(k,5)
      end do

      ! LU decomposition
      ifail = -1
      call  PNTDIAG_LU (ny,
     &     diags(1,1),diags(1,2),diags(1,3),diags(1,4),diags(1,5),
     &     ifail)


      if (ifail.ne.-1) then
         write(*,*) 'Error in pntdiag_lu. Stop.'
         stop
      end if

      if (my_node.eq.0) then
      write (*,*) '* Primary filter initialized.'
      end if

      end subroutine init_filter


c***********************************************************************


      subroutine decfac(iord,dc)
c     computes factors for (N=iord)
c
c     Qn = sum(k=0..N) (1-G)**k
c        = dc(2) + dc(3)*G + dc(4)*G**2 + ... + dc(N+2)*G**N
c        = sum(k=0..N) dc(k+2) * G**k

      implicit none
      integer, intent(in)  :: iord
      real,    intent(out) :: dc(*)
      integer              :: i,j

      do j = 1, iord + 2
         dc(j) = 0.
      end do
      dc(1) = 1
      do j = 2, iord + 2
         do i = j, 2, -1
            dc(i) = dc(i) + dc(i-1)
         end do
      end do
      dc(1) = 0
      do i = 3, iord+2, 2
         dc(i) = -dc(i)
      end do

      end subroutine decfac


c***********************************************************************


      subroutine coeff_filter_centered (dmm,dm,dp,dpp,cutoff,gew1,
     &     gew2,i3mom)
      implicit none
      real,    intent(in)  :: dmm,dm,dp,dpp,cutoff
      integer, intent(in)  :: i3mom
      real,    intent(out) :: gew1(5), gew2(5)

      real s1,s2,s3,s4
      real apx,appx,amx,ammx,anx
      real ap,app,am,amm,an
      real pi, rk, rfact

      real gim,gre
      real a,b,c

c      write(*,*) '* i3mom = ',i3mom

c---- explicit filter
c     Note: i3mom does not make a difference for aequidistant grids,
c           since the dispersion is always zero and so is the 3rd moment.
      if (i3mom.eq.0) then
         ! minimum dispersion
         s1 = -0.5
         s3 = 2.*dpp*dmm**3*dm**2+2.*dpp**2*dmm**2*dm**2-2.*dpp*dmm**4*
     &        dm-2.*dpp**2*dmm**3*dm+2.*dpp*dmm**3*dp**2+2.*dpp**2*
     &        dmm**2*dp**2+2.*dpp*dmm**4*dp+2.*dpp**2*dmm**3*dp+dp**3*
     &        dm**2*dmm-dp*dm**4*dmm+2.*dp**2*dm**2*dmm**2+dp*dm**3*
     &        dmm**2-dp**3*dm**2*dpp+dp*dm**4*dpp+2.*dp**2*dm**2*dpp**2+
     &        dp*dm**3*dpp**2+dm*dp**4*dmm-dm**3*dp**2*dmm+dm*dp**3*
     &        dmm**2-dm*dp**4*dpp+dm**3*dp**2*dpp+dm*dp**3*dpp**2
         s4 = 1./(2.*dp*dpp**4*dm+2.*dmm**3*dp**2*dm-2.*dp**2*dm**2*
     &        dmm**2-2.*dmm*dm**4*dpp+2.*dpp**2*dmm**2*dm**2-2.*dmm**2*
     &        dp**3*dpp+2.*dmm**4*dp*dm-2.*dmm*dm**3*dpp**2+2.*dmm*
     &        dp**3*dpp**2-2.*dmm**3*dm**2*dp+2.*dmm**2*dm**3*dpp+2.*
     &        dpp**3*dm**2*dp-2.*dpp**3*dp**2*dm+2.*dpp**2*dmm**2*dp**2-
     &        2.*dmm*dp**4*dpp+4.*dpp**2*dmm**4-2.*dp**2*dm**2*dpp**2+
     &        dmm**2*dp**4+dm**2*dpp**4-2.*dmm**3*dm**3+dpp**2*dp**4+4.*
     &        dpp**4*dmm**2+dmm**2*dm**4+dpp**2*dm**4+8.*dpp**3*dmm**3+
     &        2.*dpp**3*dm**3+dmm**4*dm**2+2.*dmm**3*dp**3-2.*dpp**3*
     &        dp**3+dmm**4*dp**2+dp**2*dpp**4+4.*dmm*dp**2*dpp*dm**2+2.*
     &        dmm*dp**2*dm*dpp**2-2.*dmm*dm**2*dp*dpp**2+2.*dmm**2*dp*
     &        dpp*dm**2+4.*dmm**2*dp*dm*dpp**2-2.*dmm**2*dm*dpp*dp**2)
         s2  = s3*s4
         appx= s1*s2
      else
         ! third moment zero
         appx = 0.5*dp*dm*(-dmm*dp+dp*dm-dmm**2+dmm*dm)/dpp/(-dp**2*dmm-
     &        dpp*dp**2-dp*dmm**2+dp*dpp**2+dp*dmm*dm+dp*dm*dpp+dpp**2*
     &        dmm+dmm**2*dpp+dmm**2*dm-dm*dpp**2-dmm*dm**2-dpp*dm**2)
      end if

      apx  = -0.5*(dm**2+2.*dmm*appx*dpp-dmm*dm+2.*appx*dpp**2)/(dp**2-
     &     dm**2+dmm*dp+dmm*dm)
      amx  = 0.5-apx
      ammx = (dpp*appx+dp*apx-dm*amx)/dmm
      anx  = 1.-(appx+apx+amx+ammx)

c---- implicit filter
      pi    = 4.*atan(1.)
      rk    = pi*cutoff

c     Approximation (see Diss Steffen)
c      rfact = 1./((ammx+appx)*(cos(-2.*rk)-1.)
c     &     + (amx+apx)*(cos(-rk)-1.)) + 2.

      gre = ammx*cos(-2.*rk)+amx*cos(-rk)+anx+
     &     apx*cos(rk)+appx*cos(2.*rk)
      gim = ammx*sin(-2.*rk)+amx*sin(-rk)+
     &     apx*sin(rk)+appx*sin(2.*rk)

      a = (gre-1)**2+gim**2
      b = 2*(gre-1)
      c = 1-4*gre**2-4*gim**2

      rfact = (-b-sqrt(b**2 - 4*a*c))/(2*a)


      ap    = rfact*apx
      app   = rfact*appx
      amm   = rfact*ammx
      am    = rfact*amx
      an    = 1.-(app+ap+am+amm)



c---- put coefficients together
      ! explicit filter coefficients
      gew1(1) = ammx
      gew1(2) = amx
      gew1(3) = anx
      gew1(4) = apx
      gew1(5) = appx
      ! implicit filter coefficients
      gew2(1) = amm
      gew2(2) = am
      gew2(3) = an
      gew2(4) = ap
      gew2(5) = app

      end subroutine coeff_filter_centered


c***********************************************************************


      subroutine coeff_filter_r3_l1
     &     (dmmm,dmm,dm,dp,cutoff,gew1,gew2,i3mom)

      implicit none
      integer, intent(in) :: i3mom
      real, intent(in)    :: dmmm,dmm,dm,dp,cutoff
      real, intent(out)   :: gew1(5),gew2(5)

      real   :: gew1tmp(5),gew2tmp(5)
      integer k

      call coeff_filter_r1_l3(dp,dm,dmm,dmmm,cutoff,
     &     gew1tmp,gew2tmp,i3mom)

      do k=1,5
         gew1(k) = gew1tmp(6-k)
         gew2(k) = gew2tmp(6-k)
      end do

      end subroutine coeff_filter_r3_l1


c***********************************************************************


      subroutine coeff_filter_r1_l3
     &     (dm,dp,d2,d3,cutoff,gew1,gew2,i3mom)

      implicit none
      integer, intent(in) :: i3mom
      real, intent(in)    :: dm,dp,d2,d3,cutoff
      real, intent(out)   :: gew1(5),gew2(5)

      real s1,s2,s3,s4
      real am,an,ap,a2,a3
      real amx,anx,apx,a2x,a3x
      real pi,rk,rfact

      real gim,gre
      real a,b,c

c      write(*,*) '* i3mom (r1l3) = ',i3mom

c---- explicit filter
      if (i3mom.eq.0) then
         ! minimum dispersion
         s1 = 0.5
         s3 = -dm**3*dp**2*d3+2.*d2**3*dm*dp**2-2.*d2**3*dm*d3**2+dm*dp
     &        **4*d3-2.*d2**2*dp**2*d3**2-2.*d2**3*dm**2*dp-2.*d2**4*dm*
     &        dp+2.*d2**2*dm**2*dp**2+2.*d2**3*dm**2*d3-2.*d2**2*dm**2*
     &        d3**2+2.*d2**4*dm*d3-dm*dp**3*d3**2-dp*dm**4*d3+dp**3*
     &        dm**2*d3+2.*d2**3*dp**2*d3-2.*d2**4*dp*d3+2.*d2**3*dp*
     &        d3**2-2.*dp**2*dm**2*d3**2-dp*dm**3*d3**2+d2**4*dm**2-2.*
     &        d2**3*dp**3+d2**4*dp**2+d2**2*dm**4+d2**2*dp**4+dp**2*
     &        dm**4+2.*dp**3*dm**3+dm**2*dp**4+2.*d2**3*dm**3
         s4 = 1/(-2.*dm**3*dp**2*d3-4.*d2**3*dm*d3**2+2.*dm*dp**4*d3-
     &        4.*d2**2*dp**2*d3**2-2.*d3**3*dp**2*dm+4.*d2**3*dm**2*d3-
     &        4.*d2**2*dm**2*d3**2+2.*d3**3*dm**2*dp+4.*d2**4*dm*d3-2.*
     &        dm*dp**3*d3**2-2.*dp*dm**4*d3+2.*dp**3*dm**2*d3+4.*d2**3*
     &        dp**2*d3+2.*dp*d3**4*dm-4.*d2**4*dp*d3+4.*d2**3*dp*d3**2-
     &        6.*dp**2*dm**2*d3**2-2.*dp*dm**3*d3**2+2.*d2**4*dm**2-4.*
     &        d2**3*dp**3+2.*d2**4*dp**2+2.*d2**2*dm**4+2.*d2**2*dp**4+
     &        dp**2*dm**4+2.*dp**3*dm**3+dm**2*dp**4+4.*d2**3*dm**3+4.*
     &        d2**4*d3**2-8.*d2**3*d3**3+4.*d2**2*d3**4+d3**2*dp**4-2.*
     &        d3**3*dp**3+d3**2*dm**4+2.*d3**3*dm**3+dp**2*d3**4+dm**2*
     &        d3**4)
         s2 = s3*s4
         a3x= s1*s2
      else
         ! third moment zero
         a3x=0.5*dp*dm*(-dp*dm-d2*dp+d2**2+d2*dm)/(-dm**2*dp**2-dp**2*d2
     &        *d3-dp**2*d2*dm+d3**2*dp**2+dp*d2**2*d3+dp*d2**2*dm-d3**3*
     &        dp+dp*d3*d2*dm+dm**2*d2*dp-dp*dm*d3**2-d2**2*d3**2+d3**3*
     &        d2-dm*d2**2*d3+d3**3*dm-d3*d2*dm**2+d3**2*dm**2)
      end if
      apx  = 0.5*(dm**2-2.*dm**2*a3x-2.*d2*a3x*d3+d2*dm-2.*d2*dm*a3x+
     &     2.*a3x*d3**2)/(-dp**2+dm**2+d2*dp+d2*dm)
      amx  = 0.5-apx-a3x
      a2x  = -(a3x*d3+apx*dp-amx*dm)/d2
      anx  = 1.-(a3x+a2x+apx+amx)

c---- implicit filter
      pi    = 4.*atan(1.)
      rk    = pi*cutoff

c     Approximation (see Diss Steffen)
c      ra    = a3x*cos(rk*3.)+a2x*cos(rk*2.)+anx+apx*cos(rk)+amx*cos(-rk)
c      rfact = 1./(ra-1.)+2.

      gre = amx*cos(-rk)+anx+
     &     apx*cos(rk)+a2x*cos(2.*rk)+a3x*cos(3.*rk)
      gim = amx*sin(-rk)+
     &     apx*sin(rk)+a2x*sin(2.*rk)+a3x*sin(3.*rk)

      a = (gre-1)**2+gim**2
      b = 2*(gre-1)
      c = 1-4*gre**2-4*gim**2

      rfact = (-b-sqrt(b**2 - 4*a*c))/(2*a)


      ap    = rfact*apx
      a2    = rfact*a2x
      a3    = rfact*a3x
      am    = rfact*amx
      an    = 1.-(a3+a2+ap+am)

c---- put coefficients together
      ! explicit filter coefficients
      gew1(1) = amx
      gew1(2) = anx
      gew1(3) = apx
      gew1(4) = a2x
      gew1(5) = a3x
      ! implicit filter coefficients
      gew2(1) = am
      gew2(2) = an
      gew2(3) = ap
      gew2(4) = a2
      gew2(5) = a3

      end subroutine coeff_filter_r1_l3


c***********************************************************************


      subroutine write_filter_centered (gew1,gew2,iunit)
      implicit none

      real, intent(in) :: gew1(5),gew2(5)
      integer,intent(in) :: iunit
      real ammx,amx,anx,apx,appx
      real amm,am,an,ap,app
      real pi,rk
      integer i,N
      real ra,rb,rc,rd
      real rreGexpl,rimGexpl
      real rreGimpl,rimGimpl

c      write(*,*) 'Filter weights for centered filter (expl., impl.)'
c      do k=1,5
c         write(*,*) k,gew1(k),gew2(k)
c      end do

      ammx = gew1(1)
      amx  = gew1(2)
      anx  = gew1(3)
      apx  = gew1(4)
      appx = gew1(5)

      amm  = gew2(1)
      am   = gew2(2)
      an   = gew2(3)
      ap   = gew2(4)
      app  = gew2(5)

      pi = 4.*atan(1.)
      N = 100
      do i=1,N
         rk = real(i-1)*pi/real(N-1)

c---- implicit filter
         ! Fourier transform
         ra = (ammx*cos(-rk*2.)+amx*cos(-rk)+
     &        anx+apx*cos(rk)+appx*cos(rk*2.))
         rc = (amm*cos(-rk*2.)+am*cos(-rk)+
     &        an+ap*cos(rk)+app*cos(rk*2.))
         rb = (ammx*sin(-rk*2.)+amx*sin(-rk)+
     &        apx*sin(rk)+appx*sin(rk*2.))
         rd = (amm*sin(-rk*2.)+am*sin(-rk)+
     &        ap*sin(rk)+app*sin(rk*2.))

         rreGimpl = (ra*rc+rb*rd)/(rc*rc+rd*rd)
         rimGimpl = (rb*rc-ra*rd)/(rc*rc+rd*rd)

c---- explicit filter (rc=1, rd=0)
         rc = 1.
         rd = 0.

         rreGexpl = ra
         rimGexpl = rb

         write (iunit,'(9E15.7)') rk,rreGexpl,rimGexpl,rreGimpl,rimGimpl
      end do

      end subroutine write_filter_centered


c***********************************************************************


      subroutine write_filter_r1_l3 (gew1,gew2,iunit)

      implicit none
      real, intent(in) :: gew1(5),gew2(5)
      integer, intent(in) :: iunit
      real amx,anx,apx,a2x,a3x
      real am,an,ap,a2,a3
      real pi
      real ra,rb,rc,rd,rk
      integer N,i
      real rreGexpl,rimGexpl
      real rreGimpl,rimGimpl

      amx = gew1(1)
      anx = gew1(2)
      apx = gew1(3)
      a2x = gew1(4)
      a3x = gew1(5)

      am = gew2(1)
      an = gew2(2)
      ap = gew2(3)
      a2 = gew2(4)
      a3 = gew2(5)

      pi = 4.*atan(1.)
      N = 100
      do i=1,N
         rk = (i-1)*pi/(N-1)
c---- implicit filter
         ! Fourier transform
         rc = a3*cos(rk*3.)+a2*cos(rk*2.)+
     &        an+ap*cos(rk)+am*cos(-rk)
         rd = a3*sin(rk*3.)+a2*sin(rk*2.)+
     &        ap*sin(rk)+am*sin(-rk)
         ra = a3x*cos(rk*3.)+a2x*cos(rk*2.)+
     &        anx+apx*cos(rk)+amx*cos(-rk)
         rb = a3x*sin(rk*3.)+a2x*sin(rk*2.)+
     &        apx*sin(rk)+amx*sin(-rk)

         rreGimpl = (ra*rc+rb*rd)/(rc*rc+rd*rd)
         rimGimpl = (rb*rc-ra*rd)/(rc*rc+rd*rd)

c---- explicit filter (rc=1, rd=0)
         rc = 1.
         rd = 0.

         rreGexpl = ra
         rimGexpl = rb

         write (iunit,'(9E15.7)') rk,rreGexpl,rimGexpl,rreGimpl,rimGimpl
      end do

      end subroutine write_filter_r1_l3


c***********************************************************************


      subroutine write_filter_r3_l1 (gew1,gew2,iunit)

      implicit none
      real, intent(in) :: gew1(5),gew2(5)
      integer, intent(in) :: iunit
      real ammmx,anx,ammx,amx,apx
      real ammm,an,amm,am,ap
      real pi
      real ra,rb,rc,rd,rk
      integer N,i
      real rreGexpl,rimGexpl
      real rreGimpl,rimGimpl

      ammmx = gew1(1)
      ammx = gew1(2)
      amx = gew1(3)
      anx = gew1(4)
      apx = gew1(5)

      ammm = gew2(1)
      amm = gew2(2)
      am = gew2(3)
      an = gew2(4)
      ap = gew2(5)

      pi = 4.*atan(1.)
      N = 100
      do i=1,N
         rk = (i-1)*pi/(N-1)
c---- implicit filter
         ! Fourier transform
         rc = ammm*cos(-rk*3.)+amm*cos(-rk*2.)+
     &        am*cos(-rk)+an+ap*cos(rk)
         rd = ammm*sin(-rk*3.)+amm*sin(-rk*2.)+
     &        am*sin(-rk)+ap*sin(rk)
         ra = ammmx*cos(-rk*3.)+ammx*cos(-rk*2.)+
     &        amx*cos(-rk)+anx+apx*cos(rk)
         rb = ammmx*sin(-rk*3.)+ammx*sin(-rk*2.)+
     &        amx*sin(-rk)+apx*sin(rk)


         rreGimpl = (ra*rc+rb*rd)/(rc*rc+rd*rd)
         rimGimpl = (rb*rc-ra*rd)/(rc*rc+rd*rd)

c---- explicit filter (rc=1, rd=0)
         rc = 1.
         rd = 0.

         rreGexpl = ra
         rimGexpl = rb

         write (iunit,'(9E15.7)') rk,rreGexpl,rimGexpl,rreGimpl,rimGimpl
      end do

      end subroutine write_filter_r3_l1
