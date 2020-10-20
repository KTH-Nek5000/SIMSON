c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rescale(tmax,dt,dtmax,rot,tsave,nsave,dstar,xl,
     &     fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7,
     &     fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,fpdds6,
     &     wpds1,wpds2,wpds3,wpds4,wpds5,wpds6,wpds7,wpds8,wpds9,
     &     wpds10,wpds11,wpds12,wpdds1,wpdds2,
     &     tamps,tampt,txsc,tysc,tdt,txys,tx0,
     &     fmax,fstart,fend,frise,ffall,tsmst,tsmend,tsmoo,
     &     sfd_delta,sfd_chi,dtheta0_low,dtheta0_upp,xysp,
     &     ttab,dttab)
c
c     Rescale to internal channel coordinates
c     dstar = the displacement thickness at t=0,x=0 in half channel heights
c
c     Time scale:         dstar         tcode = tphys * dstar
c     Space scale:        dstar         Lcode = Lphys * dstar
c     Velocity scale:     1             ucode = uphys
c     Acceleration scale: 1/dstar       acode = aphys / dstar
c
      implicit none

      include 'par.f'
      integer nsave
      real tmax,dt,dtmax,tsave(nsave),rot,dstar,xl
      real fpds1,fpds2,fpds3,fpds4,fpds5,fpds6,fpds7
      real fpds8,fpdds1,fpdds2,fpdds3,fpdds4,fpdds5,fpdds6
      real wpds1,wpds2,wpds3,wpds4,wpds5,wpds6,wpds7,wpds8,wpds9
      real wpds10,wpds11,wpds12
      real wpdds1,wpdds2
      real tdt,tamps,tampt,tysc,txsc,txys,tx0
      real fmax,fstart,fend,frise,ffall
      real tsmst(2),tsmend(2),tsmoo(4)
      real sfd_chi,sfd_delta
      real dtheta0_upp(scalar),dtheta0_low(scalar)
      real xysp
      real ttab,dttab

      integer i
c
c     Time stepping
c
      tmax=tmax*dstar
      dt=dt*dstar
      dtmax=dtmax*dstar
      do i=1,nsave
         tsave(i)=tsave(i)*dstar
      end do
      xysp = xysp*dstar
c
c     Rotation rate
c
      rot=rot/dstar
c
c     Localized disturbance force
c
      fpds1=fpds1*dstar
      fpds2=fpds2*dstar
      fpds3=fpds3*dstar
      fpds4=fpds4*dstar
      fpds5=fpds5*dstar
      fpds6=fpds6*dstar
      fpds7=fpds7*dstar
      fpds8=fpds8*dstar
      fpdds1=fpdds1/dstar
      fpdds2=fpdds2/dstar
      fpdds3=fpdds3/dstar
      fpdds4=fpdds4/dstar
      fpdds5=fpdds5/dstar
      fpdds6=fpdds6/dstar
c
c     Trip disturbance
c
      tamps=tamps/dstar
      tampt=tampt/dstar
      txsc=txsc*dstar
      tx0=tx0*dstar
      tysc=tysc*dstar
      tdt=tdt*dstar
c
c     Statistics sampling
c
      txys=txys*dstar
c
c     Fringe parameters
c
      fstart=fstart*dstar
      fend=fend*dstar
      frise=frise*dstar
      ffall=ffall*dstar
      fmax=fmax/dstar
c
c     Lengthening of box
c
      xl=xl*dstar
c
c     Blowing and suction
c
      wpds1=wpds1*dstar
      wpds2=wpds2*dstar
      wpds3=wpds3*dstar
      wpds4=wpds4*dstar
      wpds5=wpds5*dstar
      wpds6=wpds6*dstar
      wpds7=wpds7*dstar
      wpds8=wpds8*dstar
      wpds9=wpds9*dstar
      wpds10=wpds10*dstar
      wpds11=wpds11*dstar
      wpds12=wpds12*dstar
      wpdds1=wpdds1/dstar
      wpdds2=wpdds2/dstar
c
c     Streaks
c
      tsmst  = tsmst  * dstar
      tsmend = tsmend * dstar
      tsmoo  = tsmoo  * dstar
c
c     SFD
c
      sfd_delta = sfd_delta * dstar
      sfd_chi   = sfd_chi   / dstar

c
c     Derivative of scalar field
c
      do i=1,scalar
         dtheta0_upp(i)=dtheta0_upp(i) / dstar
         dtheta0_low(i)=dtheta0_low(i) / dstar
      end do
c
c     Changing base flow
c
      ttab=ttab*dstar
      dttab=dttab*dstar
      

      end subroutine rescale
