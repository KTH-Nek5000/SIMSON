c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cbflow(bu1,bu2,xsc,x0,eta,
     &     fbla,rlam,m1,spanv,re,nnx,xl,dstar,
     &     bstart,blength,tabfre,namfre,u0low,w0low,prex,prey,wbr,wbi,
     &     fltype,my_node,suction,asbl,vsuc,nbla,dybla,thgrad,
     &     dtheta0_low,tbc,theta0_low,theta0_upp,gr,
     &     mftab,x0tab,x0PG,ifmchange)
c
c     Construct base flow from an interpolated profile table and power-law
c     or spline interpolated freestream velocity distribution
c
c     fbla contains the Blasius solution
c     The contents of the fbla array are as follows:
c     fbla(1:nbla, 1)    f
c     fbla(1:nbla, 2)    f'
c     fbla(1:nbla, 3)    f''
c     fbla(1:nbla, 4)    f'''
c     fbla(1:nbla, 5)    g      (only if crossflow active)
c     fbla(1:nbla, 6)    g'     (only if crossflow active)
c     fbla(1:nbla, 7)    g''    (only if crossflow active)
c     fbla(1:nbla, 8)    th     (only if scalar=1)
c     fbla(1:nbla, 9)    th'    (only if scalar=1)
c     fbla(1:nbla,10)    th''   (only if scalar=1)
c
      implicit none

      include 'par.f'

      logical tabfre
      character*80 namfre
      integer fltype,nnx,nbla
      real fbla(mbla,7+3*scalar)
      real gr(scalar)
      real ybl,etabl,df1,df2,df3
      real eta(nyp),etab(nyp)
      real rlam,spanv,m1(scalar),dybla
      real x0,re,xl,dstar,xsc,bstart,blength,u0low,w0low
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real prex(nxp+15),prey(nyp*2+15)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)

      integer x,y,i,ith
      real xblu1(nxp/2),xblu2(nxp/2),xble1(nxp/2),xble2(nxp/2)
      real xblu,xble,xmin,xc
      real uu1(nxp/2),uu2(nxp/2),ue1(nxp/2),ue2(nxp/2)
      real ymu1(nxp/2),ymu2(nxp/2),yme1(nxp/2),yme2(nxp/2)
      real bstu(nxp/2),bste(nxp/2)
      real bcu(nxp/2),bce(nxp/2)
      real tmp,thuu1,thuu2,thgrad(scalar),dtheta0_low(scalar)
      real theta0_low(scalar),theta0_upp(scalar)
      integer tbc(scalar)
      real mftab,x0tab,x0PG
      logical ifmchange
c
c     Asymptotic suction layer
c
      real vsuc
      logical suction,asbl
c
c     Functions
c
      real, external :: step,xlim,cubip
c
c     Constant
c
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     MPI
c
      integer my_node
c
c     Calculate base flow in bu(x,y,1..3+scalar)
c
      if (fltype.eq.3.or.fltype.eq.9.or.fltype.eq.-1.or.fltype.eq.-2)
     &     then
c
c     Temporal boundary layers
c
         if (scalar.gt.0) call stopnow(34242)
         if (my_node.eq.0) then
            write(*,*)
            write(*,*) 'Constructing parallel base flow'
         end if
c
c     Similarity solution
c
         do y=1,nyp
            ybl = 1.+eta(y)
            etabl=ybl*sqrt(re/(2.*x0))
            df1 =       cubip(etabl,fbla(1,2),dybla,nbla)+u0low
            df2 = spanv*cubip(etabl,fbla(1,5),dybla,nbla)+w0low
            do i=1,scalar
               df3 =    cubip(etabl,fbla(1,8+3*(i-1)),dybla,nbla)
            end do
            do x=1,nnx/2
               bu1(x,y,1) = df1
               bu2(x,y,1) = df1
               bu1(x,y,2) = 0.0
               bu2(x,y,2) = 0.0
               bu1(x,y,3) = df2
               bu2(x,y,3) = df2
               do i=1,scalar
                  bu1(x,y,3+i) = df3
                  bu2(x,y,3+1) = df3
               end do
            end do
         end do
         return
      else if (fltype.eq.1.or.fltype.eq.4) then
c
c     Poiseuille base flow
c
         do y = 1,nyp
            do x = 1,nnx/2
               bu1(x,y,1) = 1. - eta(y)**2
               bu2(x,y,1) = 1. - eta(y)**2
               bu1(x,y,2) = 0.
               bu2(x,y,2) = 0.
               bu1(x,y,3) = 0.
               bu2(x,y,3) = 0.
               do ith = 1,scalar
                  bu1(x,y,3+ith) = 1.-(1.-eta(y))/2.
                  bu2(x,y,3+ith) = 1.-(1.-eta(y))/2.
               end do
            end do
         end do

         return
      else if (fltype.eq.2.or.fltype.eq.5) then
c
c     Couette base flow
c
         do y = 1,nyp
            do x = 1,nnx/2
               bu1(x,y,1) = eta(y)
               bu2(x,y,1) = eta(y)
               bu1(x,y,2) = 0.
               bu2(x,y,2) = 0.
               bu1(x,y,3) = 0.
               bu2(x,y,3) = 0.
               do ith = 1,scalar
                  bu1(x,y,3+ith) = 1.-(1.-eta(y))/2.
                  bu2(x,y,3+ith) = 1.-(1.-eta(y))/2.
               end do
            end do
         end do
         return
      else if (abs(fltype).eq.20) then
c
c     Buoyancy-driven boundary layer flow
c
         do y=1,nyp
            etab(y) = (eta(y)+1.)/dstar
         end do
         do y = 1,nyp
            do x = 1,nnx/2
               bu1(x,y,1) = gr(1)/2*exp(-etab(y))*sin(etab(y))
               bu2(x,y,1) = gr(1)/2*exp(-etab(y))*sin(etab(y))
               bu1(x,y,2) = 0.
               bu2(x,y,2) = 0.
               bu1(x,y,3) = 0.
               bu2(x,y,3) = 0.
               bu1(x,y,4) = exp(-etab(y))*cos(etab(y))
               bu2(x,y,4) = exp(-etab(y))*cos(etab(y))
            end do
         end do
         return
      else if (fltype.eq.-3) then
c
c     Asymptotic suction boundary layer
c
         if (.not.asbl) then
            if (my_node.eq.0) then
               write(*,*)
     &              'Base flow not known for this suction rate: ',
     &              vsuc
            end if
            call stopnow(34324)
         end if

         if (scalar.ge.1) then
            if (my_node.eq.0) then
               write(*,*) 'ASBL with scalar not implemented.'
            end if
            call stopnow(5435)
         end if

         do y=1,nyp
            df1=1.-exp(-(1.+eta(y))/dstar)+u0low
            do x=1,nnx/2
               bu1(x,y,1)=df1
               bu2(x,y,1)=df1
               bu1(x,y,2)=-vsuc
               bu2(x,y,2)=-vsuc
               bu1(x,y,3)=0.
               bu2(x,y,3)=0.
            end do
         end do

         return

      end if
c
c     Here we continue for fltype=6,7,8 (Blasius/Falkner-Skan-Cooke)
c
      if (my_node.eq.0) then
         write(*,*)'Construct non-parallel base flow with m=',rlam
      end if
      do x=1,nnx/2
c
c     By blending the flow from the upstream and the downstream part
c     we need a coordinate which starts at bstart and goes on
c     to fend+xl
c
c     To avoid complex looping we first construct the x-coordinates
c     for each point since we need to blend the 'upstream' and
c     'downstream'
c
c     Boundary layer we need corresponding coordinates
c     thus xblu1,xble1 are in the interval [x0+bstart,x0+bstart+xl]
c     and  xblu2,xble2 are in the interval [x0+bstart+xl,x0+bstart+2*xl]
c
c
c     Create coordinates and blending function
c     XXXu are the odd points (2*x-1)
c     XXXe are the even points (2*x)
c     XXX1 are the coordinates in interval [bstart,bstart+xl[
c     XXX2 are the coordinates in interval [bstart+xl,bstart+2*xl[
c
c     Odd points (2*x-1)
c
         xc   = real(2*x-1-nnx/2-1)/real(nnx)*xl+xsc
         xblu = xc-(int((xc-bstart)/xl+1.)-1.)*xl
c
c     bstu is the blending step function
c
         bstu(x)  = step((xblu-bstart)/blength)
         xblu1(x) = xblu+x0
         xblu2(x) = xblu1(x)+xl
c
c     Even points (2*x)
c
         xc   = real(2*x-nnx/2-1)/real(nnx)*xl+xsc
         xble = xc-(int((xc-bstart)/xl+1.)-1.)*xl
c
c     bste is the blending step function
c
         bste(x)  = step((xble-bstart)/blength)
         xble1(x) = xble+x0
         xble2(x) = xble1(x)+xl
      end do
c
c     Now limit the coordinate to avoid going upstream of leading edge
c
      xmin = x0/2.
c     changed from 0.5
c     changed from 0.72

c
c     Fill the coordinate and free-stream arrays
c
      do x=1,nnx/2
         xblu1(x)=xlim(xblu1(x),xmin)
         uu1(x)=(xblu1(x)/x0)**rlam
         ymu1(x)=sqrt((rlam+1.)*uu1(x)*re/(2.*xblu1(x)))

         xble1(x)=xlim(xble1(x),xmin)
         ue1(x)=(xble1(x)/x0)**rlam
         yme1(x)=sqrt((rlam+1.)*ue1(x)*re/(2.*xble1(x)))

         xblu2(x)=xlim(xblu2(x),xmin)
         uu2(x)=(xblu2(x)/x0)**rlam
         ymu2(x)=sqrt((rlam+1.)*uu2(x)*re/(2.*xblu2(x)))

         xble2(x)=xlim(xble2(x),xmin)
         ue2(x)=(xble2(x)/x0)**rlam
         yme2(x)=sqrt((rlam+1.)*ue2(x)*re/(2.*xble2(x)))
      end do
c
c     For spline interpolated streamwise freestream velocity
c     read the velocity table and interpolate onto uu and ue
c
      if (tabfre) then
         call splfre(uu1,ue1,uu2,ue2,bstart,dstar,xl,xsc,namfre,
     &        mftab,x0tab,x0PG,ifmchange)
      end if
c
c     Now compute u,w,th from the similarity solution
c
      do y=1,nyp
         ybl=1.+eta(y)
         do x=1,nnx/2
c     Schlichting scaled eta from xbl1 and interpolate u
            etabl=ybl*ymu1(x)
            df1=cubip(etabl,fbla(1,2),dybla,nbla)
            etabl=ybl*ymu2(x)
            df2=cubip(etabl,fbla(1,2),dybla,nbla)
c     Compute streamwise base flow
            bu1(x,y,1)=uu1(x)*df1*bstu(x)+uu2(x)*df2*(1.-bstu(x))+u0low

c     Schlichting scaled eta from xbl1 and interpolate u
            etabl=ybl*yme1(x)
            df1=cubip(etabl,fbla(1,2),dybla,nbla)
            etabl=ybl*yme2(x)
            df2=cubip(etabl,fbla(1,2),dybla,nbla)
c     Compute streamwise base flow
            bu2(x,y,1)=ue1(x)*df1*bste(x)+ue2(x)*df2*(1.-bste(x))+u0low

c     Schlichting scaled eta from xbl1 and interpolate u
            etabl=ybl*ymu1(x)
            df1=cubip(etabl,fbla(1,5),dybla,nbla)
            etabl=ybl*ymu2(x)
            df2=cubip(etabl,fbla(1,5),dybla,nbla)
c     Compute spanwise base flow
            bu1(x,y,3)=spanv*(df1*bstu(x)+df2*(1.-bstu(x)))+w0low

c     Schlichting scaled eta from xbl1 and interpolate u
            etabl=ybl*yme1(x)
            df1=cubip(etabl,fbla(1,5),dybla,nbla)
            etabl=ybl*yme2(x)
            df2=cubip(etabl,fbla(1,5),dybla,nbla)
c     Compute spanwise base flow
            bu2(x,y,3)=spanv*(df1*bste(x)+df2*(1.-bste(x)))+w0low

            do i=1,scalar
c     Schlichting scaled eta from xbl1 and interpolate u
               etabl=ybl*ymu1(x)
               thuu1=(xblu1(x)/x0)**m1(i)
               df1=cubip(etabl,fbla(1,8+3*(i-1)),dybla,nbla)
               etabl=ybl*ymu2(x)
               thuu2=(xblu2(x)/x0)**m1(i)
               df2=cubip(etabl,fbla(1,8+3*(i-1)),dybla,nbla)
c     Compute scalar base flow
               bu1(x,y,3+i)=thuu1*df1*bstu(x)+thuu2*df2*(1.-bstu(x))

c     Schlichting scaled eta from xbl1 and interpolate u
               etabl=ybl*yme1(x)
               thuu1=(xble1(x)/x0)**m1(i)
               df1=cubip(etabl,fbla(1,8+3*(i-1)),dybla,nbla)
               etabl=ybl*yme2(x)
               thuu2=(xble2(x)/x0)**m1(i)
               df2=cubip(etabl,fbla(1,8+3*(i-1)),dybla,nbla)
c     Compute scalar base flow
               bu2(x,y,3+i)=thuu1*df1*bste(x)+thuu2*df2*(1.-bste(x))
            end do
         end do
      end do
c
c     Now construct the v velocity from u by continuity
c
      do y=1,nyp
         do x=1,nnx/2
            bu1(x,y,2)=bu1(x,y,1)
            bu2(x,y,2)=bu2(x,y,1)
         end do
      end do
c
c     First Fourier transform
c
      call vrfftf(bu1(1,1,2),bu2(1,1,2),wbr,wbi,nnx,nyp,1,nxp/2+1,prex)
c
c     Set odd-ball to zero
c
      do y=1,nyp
         bu1(nnx/2+1,y,2)=0.0
      end do
c
c     Chebyshev transform
c
      call vchbf(bu1(1,1,2),wbr,nyp,nnx/2,nxp/2+1,1,prey)
      call vchbf(bu2(1,1,2),wbr,nyp,nnx/2,nxp/2+1,1,prey)
c
c     Normalize
c
      do y=1,nyp
         do x=1,nnx/2
            bu1(x,y,2)=bu1(x,y,2)*(2./real(nnx*(nyp-1)))
            bu2(x,y,2)=bu2(x,y,2)*(2./real(nnx*(nyp-1)))
         end do
      end do
c
c     Compute dvdy = -dudx = -ialpha u
c
      do y=1,nyp
         do x=1,nnx/2
            tmp       =  2.*pi/xl*real(x-1)*bu2(x,y,2)
            bu2(x,y,2)= -2.*pi/xl*real(x-1)*bu1(x,y,2)
            bu1(x,y,2)= tmp
         end do
      end do
c
c     Then integrate
c
      do x=1,nnx/2
         bcu(x)=0.0
         bce(x)=0.0
      end do
      call icheb(wbr,bu1(1,1,2),bcu,nyp,nnx/2,nxp/2+1)
      call icheb(wbi,bu2(1,1,2),bce,nyp,nnx/2,nxp/2+1)
      do y=1,nyp
         do x=1,nnx/2
            bu1(x,y,2)=wbr(x,y)
            bu2(x,y,2)=wbi(x,y)
         end do
      end do
c
c     Now find the value at the lower boundary for each mode
c     note skipping of mode 1
c
      do y=2,nyp-1,2
         do x=1,nnx/2
            bcu(x) = bcu(x)+bu1(x,y,2)-bu1(x,y+1,2)
            bce(x) = bce(x)+bu2(x,y,2)-bu2(x,y+1,2)
         end do
      end do
c
c     To set bc v=0 for y=0
c
      do x=1,nnx/2
         bu1(x,1,2) = bcu(x)
         bu2(x,1,2) = bce(x)
      end do
c
c     Go back to physical space
c     Chebyshev transform
c
      call vchbb(bu1(1,1,2),wbr,nyp,nnx/2,nxp/2+1,1,prey)
      call vchbb(bu2(1,1,2),wbr,nyp,nnx/2,nxp/2+1,1,prey)
c
c     Fourier transform
c
      call vrfftb(bu1(1,1,2),bu2(1,1,2),wbr,wbi,nnx,nyp,1,nxp/2+1,prex)
c
c     Rescale scalar base flow depending on boundary condition
c
      do i=1,scalar
         if (my_node.eq.0) then
            write(*,*) 'Scalar ',i,' m1=',m1(i),' tbc = ',tbc(i)
         end if
         if (abs(m1(i)-0.5).lt.1e-13) then
            if (my_node.eq.0) then
               write(*,*) 'm1=0.5 detected --> constant wall-gradient'
               write(*,*) '  thgrad (old) = ',thgrad(i)
               write(*,*) '  thgrad (new) = ',dtheta0_low(i)*dstar
               write(*,*) 'freestream val.= ',theta0_upp(i)
            end if
            do y=1,nyp
               do x=1,nnx/2
                  bu1(x,y,3+i) = theta0_upp(i)+bu1(x,y,3+i) *
     &                 (dtheta0_low(i)/thgrad(i)*dstar)
                  bu2(x,y,3+i) = theta0_upp(i)+bu2(x,y,3+i) *
     &                 (dtheta0_low(i)/thgrad(i)*dstar)
               end do
            end do
         end if
         if (m1(i).eq.0) then
            if (my_node.eq.0) then
               write(*,*) 'm1=0 detected'
               write(*,*) 'wall value       : ',theta0_low(i)
               write(*,*) 'freestream value : ',theta0_upp(i)
            end if
            do y=1,nyp
               do x=1,nnx/2
                  bu1(x,y,3+i) = theta0_upp(i) +
     &                 (theta0_low(i)-theta0_upp(i))*bu1(x,y,3+i)
                  bu2(x,y,3+i) = theta0_upp(i) +
     &                 (theta0_low(i)-theta0_upp(i))*bu2(x,y,3+i)
               end do
            end do
         end if
      end do

      end subroutine cbflow
