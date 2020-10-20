c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine bflow(bu1,bu2,bu1jr,bu1ji,xs0,x0,eta,
     &     fbla,rlam,m1,spanv,re,xl,dstar,
     &     bstart,blength,tabfre,namfre,rbfl,nambfl,
     &     u0low,w0low,u0upp,du0upp,wd1,prex,prey,
     &     wbr,wbi,boxr,boxi,w3,fltype,my_node,
     &     suction,asbl,vsuc,pert,nbla,dybla,thgrad,dtheta0_low,tbc,
     &     theta0_low,theta0_upp,w0upp,bf3,gr,
     &     mftab,x0tab,x0PG,ifmchange)
c
c     Get base flow to be used for fringe forcing, temporal forcing and
c     and compute normal derivatives of the base flow at the upper
c     boundary
c
      implicit none

      include 'par.f'

      real gr(scalar)
      logical tabfre,rbfl
      character*80 namfre,nambfl
      real fbla(mbla,7+3*scalar)
      integer fltype
      real rlam,m1(scalar),spanv,w0upp
      real x0,re,xl,dstar,xs0,bstart,blength,u0low,w0low,du0upp,u0upp
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bu1jr(nxp/2+1,3+scalar,3),bu1ji(nxp/2+1,3+scalar,3)
      real wd1(nyp,4)
      real prex(nxp+15),prey(nyp*2+15)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real w3(nx/2,mbz,nyp)
      logical suction,asbl
      real vsuc
      real eta(nyp)
      real nbla,dybla,thgrad(scalar),dtheta0_low(scalar)
      real theta0_low(scalar),theta0_upp(scalar)
      integer tbc(scalar)
      real mftab,x0tab,x0PG

      integer x,y,i,j,my_node

      logical pert, bf3,ifmchange

c     in order to force reading a base flow file
c      rbfl = .true.
c      nambfl = 'ts520.u'

      if (rbfl.and..not.bf3) then
c
c     Read the base flow from a standard 3D flow field
c
         if (my_node.eq.0) then
            write(ios,*) 'Read basic flow from 3D field ',trim(nambfl)
         end if

         call rbflow(bu1,bu2,nambfl,u0low,w0low,prex,wbr,wbi,
     &        boxr,boxi,w3)
      else if (.not.bf3) then
c
c     Compute base flow from boundary layer profiles and power law
c     or spline interpolated velocity distribution at upper wall
c
         call cbflow(bu1,bu2,xs0,x0,eta,
     &        fbla,rlam,m1,spanv,re,nxp,xl,dstar,
     &        bstart,blength,tabfre,namfre,u0low,w0low,prex,prey,
     &        wbr,wbi,fltype,my_node,suction,asbl,vsuc,nbla,dybla,
     &        thgrad,dtheta0_low,tbc,theta0_low,theta0_upp,gr,
     &        mftab,x0tab,x0PG,ifmchange)
      end if
c
c     Now compute the value and derivative of all functions along the upper
c     boundary:
c     bu1jr(x,u i=1/v i=2/w i=3,u/v/w j=1,dudy etc j=2,d2udy2 etc j=3)
c
      do j=1,3
         do i=1,3+scalar
            do x=1,nxp/2+1
               bu1jr(x,i,j)=0.0
               bu1ji(x,i,j)=0.0
            end do
         end do
      end do
      do i=1,3+scalar
         do x=1,nxp/2
            bu1jr(x,i,1)=bu1(x,1,i)
            bu1ji(x,i,1)=bu2(x,1,i)
         end do
      end do
      do i=1,3+scalar
         do y=1,nyp
            do x=1,nxp/2
               bu1jr(x,i,2)=bu1jr(x,i,2)+bu1(x,y,i)*wd1(y,1)
               bu1ji(x,i,2)=bu1ji(x,i,2)+bu2(x,y,i)*wd1(y,1)
               bu1jr(x,i,3)=bu1jr(x,i,3)+bu1(x,y,i)*wd1(y,2)
               bu1ji(x,i,3)=bu1ji(x,i,3)+bu2(x,y,i)*wd1(y,2)
            end do
         end do
      end do
      do j=1,3
         do i=1,3+scalar
c
c     Fourier transform derivatives and normalize
c
            call vrfftf(bu1jr(1,i,j),bu1ji(1,i,j),
     &           wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               bu1jr(x,i,j)=bu1jr(x,i,j)*(1./real(nxp))
               bu1ji(x,i,j)=bu1ji(x,i,j)*(1./real(nxp))
            end do
         end do
      end do

      if (pert) then
         do j=1,3
            do i=1,3+scalar
               do x=1,nxp/2+1
                  bu1jr(x,i,j)=0.0
                  bu1ji(x,i,j)=0.0
               end do
            end do
         end do
         u0low = 0.
         w0low = 0.
         w0upp = 0.
         u0upp = 0.
      end if
c
c     Add periodic point in physical space
c
      do i=1,3+scalar
         do y=1,nyp
            bu1(nxp/2+1,y,i)=bu1(1,y,i)
            bu2(nxp/2+1,y,i)=bu2(1,y,i)
         end do
      end do
c
c     Set the upper values according to base flow at x=0 (inlet)
c
      u0upp  = bu1jr(1,1,1)
      du0upp = bu1jr(1,1,2)
      w0upp  = bu1jr(1,3,1)


      end subroutine bflow
