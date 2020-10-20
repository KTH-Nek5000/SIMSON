c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program pxyst
c
c     xy statistixs plotting
c
c
      implicit none

      include 'par.f'

      integer mxys,nxys2,mgr,nxysth2,mxysth,scalarind
c
c     When changing nxys in par.f,
c     or nxys2 or nxysth2 (number of derived statistics)
c     note that the number of data elements in the
c     data statement for pdim,pdim2,thdim and thdim2 must be changed below
c
      parameter (nxys2=274)
      parameter (nxysth2=100)
      parameter (mgr=5)
      real totxys(nx,nyp,nxys),pxys(nx,nyp,nxys)
      real pxys2(nx,nyp,nxys2)
      real totxysth(nx,nyp,nxysth,scalar),pxysth(nx,nyp,nxysth)
      real pxysth2(nx,nyp,nxysth2)
      real utau(nx),uinf(nx),delta(nx,7)
      real nux(nx),tinf(nx),twall(nx),thdelta(nx,2)
      real deltar(nx),beta(nx)
      real grid(nx+2,nyp,2),plane(nx+2,nyp),w(nx,nyp)
      real xgraph(nx+1,mgr),ygraph(nyp,mgr)
      real xgrid(nx+1,mgr),ygrid(nyp,mgr)
      integer pdim(nxys,2),pdim2(nxys2,2)
      integer thdim(nxysth,3),thdim2(nxysth2,3)

      equivalence(totxys,pxys)

      character*1 ans
      character*80 namxys,xmat
      real mhd_n
c
c     Geometrics
c
      real eta(nyp)
      real pi,xl,zl,dx,xs,xr,re,t,dstar,xrd,pr(scalar)
      parameter (pi = 3.1415926535897932385)
c
c     Flow
c
      integer fltype
      real blength,bstart,rlam,spanv,sumw
c
c     Plot parameters
c
      real xp(mgr),yp(mgr)
      integer x1(mgr),y1(mgr),npy,npx,i
      logical lpgmr,ltab,lmat
      integer mx,my,ivar,iplot,iscale,ifilt,ivarth
      integer mbox,ncont,x,y,iopt
      integer mxr,npoint(mgr)
      real blev,elev,dlev,sfilt
      character*80 heada,headb,headc
      character*7 var,xscale
      character*12 scl
      integer ihard,firstmax
c
c     Initialize physical dimensions of each variable :
c     each variable in pxys or pxys2 is described as
c     (velocity)^pdim[2](j,1)*(length)^pdim[2](j,2)
c     note order : first velocity dimension of all variables
c     then length dimensions of all variables
c     each variable in pxysth or pxysth2 is described as
c     (velocity)^thdim[2](j,1)*(length)^thdim[2](j,2)*(scalar)^thdim[2](j,3)
c     note order : first velocity dimension of all variables
c     then length dimensions of all variables
c     lastly scalar dimensions of all variables
c
      data pdim/
c 001-052
     &     6*1,6*1,3*2,9*2,6*2 ,2*2,3*3,7*3 ,2*1, 2,6*1 ,0,
c 053-082
     &     3 ,6*2,1,2*3 ,10*3,2*3 ,3*4,6,8,3*3,
c 083-096
     &     2*1,6*1,3*2,3*2,

     &     6*0,6*-1,3*0,9*0,6*-2,2*0,3*0,7*-1,2*0,-2,6*-1,0,
     &     -1,6*0,1,2*-1,10*0,2*-1,3*0,0,0,3*0,
     &     2*1,6*0,3*0,3*1/
      data pdim2/
c 001-044
     &     2,3*2, 3*2, 0,0,4*3, 6*1, 7*3 ,0,0,2 ,3 ,2 ,3 ,10*3 ,2 ,4 ,
c 045-085
     &     2 ,7*2 ,6*2 ,2*2 ,6*1 ,2*0,2*0,1,0,0,0,1,1,0,0,1,0,0,3*3 ,0,
c 086-141
     &     3*3 ,3 ,0,3*1,2,3*2,2,2 ,2*2 ,4*2 ,4*2 ,5*3 ,2*2 ,8*0,17*3,
c 142-268
     &     123*3 ,3*3,2,
c     269-272
     &     6*3,
c
     &     0,3*-2,3*-2,1,0,4*-1,6*-1,7*-1,0,1,-1,-2,-1,-2,10*-1,-1,-2,
     &     -3,7*-2,6*-1,2*-2,6*-2,2*1,2*0,0,1,0,0,0,0,0,0,0,0,1,3*-1,0,
     &     3*-1,-1,0,3*0,0,3*0,0,-2,2*-1,4*-2,4*-1,5*-1,2*-1,8*0,17*-1,
     &     123*-1,3*-1,0,
     &     6*-1/



      data thdim/
c 001-035
     &     2*0,3*1,3*1,3*0 ,2,3*2 ,6*2,9*1 ,3*1 ,0,0,3*1 ,3*1,3*1,
     &     2*0,3*0,3*0,3*-2,0,3*-1,6*0,9*-1,3*-2,0,0,3*-1,3*0,3*0,
     &     2*1,3*1,3*2,3*2 ,1,3*1 ,6*1,9*1 ,3*1 ,3,4,3*2 ,3*2,3*1/
      data thdim2/
c 001-094
     &     2*0 ,6*1 ,9*1 ,3*0,0,2*0 ,2*0 ,12*1 ,57*2 ,6*1 ,
     &     2*-1,6*-1,9*-1,3*0,0,2*-1,2*-2,12*-1,57*-1,6*-1,
     &     2*1 ,6*1 ,9*1 ,3*0,2,2*2 ,2*2 ,12*2 ,57*1 ,6*2/

c-----------------------------------------------------------------------

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                        pxyst $Rev$  '//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)
c
c     Write out scalings
c
      write(11,*) '#,velocity,length'
      do i=1,nxys
         write(11,*) i,pdim(i,1),pdim(i,2)
      end do
      write(11,*)
      do i=1,nxys2
         write(11,*) -i,pdim2(i,1),pdim2(i,2)
      end do
      write(11,*) '#,velocity,length,scalar'
      do i=1,nxysth
         write(11,*) i,thdim(i,1),thdim(i,2),thdim(i,3)
      end do
      write(11,*)
      do i=1,nxysth2
         write(11,*) -i,thdim2(i,1),thdim2(i,2),thdim2(i,3)
      end do
      close(11)



      write(*,*) 'Give file to plot from'
      read(*,'(a)') namxys
      call  rxys(totxys,totxysth,mxys,mxysth,sumw,re,xl,zl,t,dstar,
     &     fltype,bstart,blength,rlam,spanv,namxys,w,pdim,thdim,pr,
     &     scalarind,mhd_n)

c
c     When reading correlation file
c
      if (re.lt.0.) then
         call correlation(namxys)
      end if
      xs=0.0
      write(*,*) '        read ',mxys,' statistics
     &for velocity and pressure'
      write(*,*) 'compiled for ',nxys,' statistics
     &for velocity and pressure'
      write(*,*) '        read ',mxysth,' statistics for scalar'
      write(*,*) 'compiled for ',nxysth,' statistics for scalar'
      write(*,*) '* time ',t
      write(*,*) '* sumw ',sumw/dstar
      do y=1,nyp
         eta(y)=cos(pi*real(y-1)/real(nyp-1))
      end do

      if (fltype.lt.0.or.fltype.eq.3.or.fltype.ge.6.and.fltype.le.8)
     &     then
         do y=1,nyp
            eta(y)=(1.+eta(y))/dstar
         end do
      end if
      write(*,9010)xl,zl,t,re
 9010 format(' xl= ',F8.2,' zl= ',f8.2,' t= ',F8.1,' re= ',F8.1)
      if (fltype.eq.3) write(*,9015) 2./dstar
 9015 format(' yl= ',F8.2)
      if (fltype.eq.-2) write(*,*) 'Falkner-Skan-Cooke boundary layer'
      if (fltype.eq.-1) write(*,*) 'Falkner-Skan boundary layer'
      if (fltype.eq.1)  write(*,*) 'Poiseuille flow'
      if (fltype.eq.2)  write(*,*) 'Couette flow'
      if (fltype.eq.3)  write(*,*) 'Blasius boundary layer flow'
      if (fltype.eq.4)  write(*,*) 'spatial Poiseuille flow'
      if (fltype.eq.5)  write(*,*) 'spatial Couette flow'
      if (fltype.eq.6)  write(*,*) 'spatial Blasius boundary layer'
      if (fltype.eq.7)  write(*,*) 'spatial Falkner-Skan boundary layer'
      if (fltype.eq.8)
     &     write(*,*) 'spatial Falkner-Skan-Cooke boundary layer'
      if (fltype.eq.9) write(*,*) 'spatial parallel boundary layer'

      write(*,*) 'filter input data (y/n)'
      read(*,'(a)') ans
      if (ans.ne.'N'.and.ans.ne.'n') then
         write(*,*) 'give filter lengthscale l'
         write(*,*) 'the filter is exp(-(k*l/2pi)**2)'
         read(*,*) sfilt
         call infilt(totxys,totxysth,sfilt,xl,scalarind)
      end if
c
c     Compute derived field variables
c
      write(*,*) 'doing pvar'
      call pvar(pxys,pxys2,pxysth,pxysth2,nxys2,nxysth2,
     &     totxys,totxysth,xl,zl,eta,re,dstar,pr,scalarind,mhd_n,
     &     firstmax)
c
c     Find local skin friction and displacement thickness
c
      write(*,*) 'doing locsc'
      call locsc(utau,uinf,delta,deltar,beta,pxys,pxys2,pxysth,
     &     re,dstar,xl,nux,thdelta,pr,scalarind,fltype,tinf,twall)
      write(*,*) 'generate mean flow plots ?'
      read(*,'(a)') ans
      if (ans.ne.'N'.and.ans.ne.'n')
     &     call meanpl(utau,uinf,delta,deltar,beta,xl,zl,re,t,namxys,
     &     nux,thdelta,pr,scalarind,fltype)
c
c     Plot loop
c
 1000 continue

      ivar = 0
      ivarth = 0

      write(*,*) 'give statistic to plot'
      write(*,*) '1-3   u,v,w'
      write(*,*) '4-6   urms,vrms,wrms **'
      write(*,*) '7-9   omx,omy,omz'
      write(*,*) '10-12 omxrms,omyrms,omzrms **'
      write(*,*) '13-15 uv,uw,vw **'
      write(*,*) '16-18 u,v,w products with x+1'
      write(*,*) '19-21 u,v,w products with y+1'
      write(*,*) '22-24 u,v,w products with z+1'
      write(*,*) '25-30 re*(eps_11,eps_22,eps_33,eps_12,eps_13,eps_23)'
      write(*,*) '31    pressure'
      write(*,*) '32    pressure_rms **'
      write(*,*) '33-35 pu,pv,pw'
      write(*,*) '36-42 p ux,p vy,p wz,p uy,p vz,p wx,p uz'
      write(*,*) '43-44 min/max u'
      write(*,*) '45    <S_ij S_ij>'
      write(*,*) '46-51 <S_ij>'
      write(*,*) '52    Smagorinsky coefficient C=C_S^2'
      write(*,*) '53    <tau_ij S_ij>'
      write(*,*) '54-59 <tau_ij>'
      write(*,*) '60    <nu_t>'
      write(*,*) '61/62 forward/backscatter'
      write(*,*) '63-65 u^3, v^3, w^3'
      write(*,*) '66-72 v u^2,w u^2,u v^2,w v^2,u w^2,v w^2,uvw'
      write(*,*) '73-74 p vx,p wy'
      write(*,*) '75-77 u^4, v^4, w^4'
      write(*,*) '78-79 p^3,p^4'
      write(*,*) '80-82 u_j tau_ij'
      write(*,*) '83-84 phi,phi_rms'
      write(*,*) '85-93 <j_k>, j_k,rms, <j1j2>,<j1j3>,<j2j3>'
      write(*,*) '94-96 <phi j_k>'
      write(*,*)
      write(*,*) "-1       turbulent kinetic energy k=1/2*<ui'ui'>"
      write(*,*) '-(2-4)   q2ro_11,q2ro_22,q2ro_33'
      write(*,*) '-(5-7)   lapl urms2,vrms2,wrms2'
      write(*,*) '-8       taylor microscale lambda_g from two-pt corr'
      write(*,*) '-9       re_lambda from two-point corr'
      write(*,*) '-(10-12) epsilon_11,22,33  from two-point corr'
      write(*,*) '-13      dissipation eps from two-point corr'
      write(*,*) '-(14-19) ux,vx,wx,uy,vy,wy'
      write(*,*) '-(20-25) eps_11,eps_22,eps_33,eps_12,eps_13,eps_23'
      write(*,*) '-26      dissipation eps'
      write(*,*) '-27      -uv/(2k)'
      write(*,*) '-28      taylor microscale lambda'
      write(*,*) '-(29-32) kx,epsx,ky,epsy'
      write(*,*) '-(33-38) (pu)x,(pv)x,(pw)x,(pu)y,(pv)y,(pw)y'
      write(*,*) '-(39-42) u px,v px,u py,v py'
      write(*,*) '-43      (urms2-vrms2)x,'
      write(*,*) '-44      eps/k*(urms2-vrms2)*ux,'
      write(*,*) '-45      (uv)yyy'
      write(*,*) '-46      (wrms2)xx'
      write(*,*) '-(47-48) (urms2)xx,(urms2)yy'
      write(*,*) '-(49-50) (vrms2)xx,(vrms2)yy'
      write(*,*) '-(51-52) (uv)xx,(uv)yy'
      write(*,*) '-(53-56) (urms2)x,(vrms2)x,(urms2)y,(vrms2)y'
      write(*,*) '-(57-58) (uv)x,(uv)y'
      write(*,*) '-(59-60) kxx,kyy'
      write(*,*) '-(61-66) uxx,vxx,wxx,uyy,vyy,wyy'
      write(*,*) '-(67-70) Lk, Lu, Lu/Lk, Lu/Lk (norm)'
      write(*,*) '-(71-72) max(urms), y(max(urms))'
      write(*,*) '-(710-720) max(urms), y(max(urms))'
      write(*,*) '-(73-74) Tu (lokal), Tu (freestream)'
      write(*,*) '-(75)    Blasius profile'
      write(*,*) '-(76-78) dev u-ublas, norm max, norm min'
      write(*,*) '-(79-80) max(rms)-min(rms), norm'
      write(*,*) '-81      delta'
      write(*,*) '-82      eps_visc,tot: -2/Re<S_ij S_ij>'
      write(*,*) '-83      eps_visc,mean: -2/Re<S_ij> <S_ij>'
      write(*,*) "-84      eps_visc,fluct: -2/Re<S_ij' S_ij'>"
      write(*,*) '-85      Smagorinsky coefficient C_S'
      write(*,*) '-86      eps_sgs,tot: <tau_ij S_ij>'
      write(*,*) '-87      eps_sgs,mean: <tau_ij> <S_ij>'
      write(*,*) "-88      eps_sgs,fluct: <tau_ij' S_ij'>"
      write(*,*) "-89      production -<u_i'u_j'><S_ij>"
      write(*,*) '-90      P/eps_fluct'
      write(*,*) "-91      sqrt(<u1'u1'>+<tau_11>)"
      write(*,*) "-92      sqrt(<u2'u2'>+<tau_22>)"
      write(*,*) "-93      sqrt(<u3'u3'>+<tau_33>)"
      write(*,*) "-94      <u1'u2'>+<tau_12>"
      write(*,*) "-95      <u1'u1'>+<tau_11>-1/3*(<ui'ui'>)"
      write(*,*) "-96      <u2'u2'>+<tau_22>-1/3*(<ui'ui'>)"
      write(*,*) "-97      <u3'u3'>+<tau_33>-1/3*(<ui'ui'>)"
      write(*,*) "-98      trace of tau_ij: <tau_kk>"
      write(*,*) '-(99-101) (wrms2)yy,(wrms2)x,(wrms2)y'
      write(*,*) '-(102-105) <uw>xx,<uw>yy,<vw>xx,<vw>yy'
      write(*,*) '-(106-109) <uw>x,<uw>y,<vw>x,<vw>y'
      write(*,*) '-(110-114) w px,w py,u pz,v pz,w pz'
      write(*,*) '-(115-116) px,py'
      write(*,*) '-(117-119) skewness for the velocities'
      write(*,*) '-(120-122) flatness for the velocities'
      write(*,*) '-(123-124) skewness and flatness for the pressure'
      write(*,*) '-125     con,tot: 1/2d<uiuiuj>/dxj'
      write(*,*) '-126     con,mean: 1/2d(<ui><ui><uj>)/dxj'
      write(*,*) "-127     con,fluct: 1/2<uj>d<ui'ui'>)/dxj"
      write(*,*) '-128     mol_dif,tot: 1/(2Re)*<d2ui^2/dxjdxj>'
      write(*,*) '-129     mol_dif,mean:1/(2Re)*d2<ui>^2/dxjdxj'
      write(*,*) "-130     mol_dif,fluct: 1/(2Re)*<d2ui^2'/dxjdxj>"
      write(*,*) "-131     tur_dif,mean:-d(<ui'uj'><ui>)/dxj"
      write(*,*) "-132     tur_dif,fluct:-1/2*d(<ui'uj'ui'>)/dxj"
      write(*,*) '-133     v-p_gra,tot: -<ui dp/dxi>'
      write(*,*) '-134     v-p_gra,mean: -<ui><dp/dxi>'
      write(*,*) "-135     v-p_gra,fluct: -<ui' dp'/dxi>"
      write(*,*) '-136     dis_eps,tot:-1/Re<dui/dxj dui/dxj>'
      write(*,*) '-137     dis_eps,mean:-1/Re<dui/dxj><dui/dxj>'
      write(*,*) "-138     dis_eps,fluct:-1/Re<dui'/dxj dui'/dxj>"
      write(*,*) '-139     vis_dif tot:2/Red<uiSij>/dxj '
      write(*,*) '-140     vis_dif mean:2/Red<ui><Sij>/dxj '
      write(*,*) "-141     vis_dif fluct:2/Red<ui'Sij'>/dxj"
      write(*,*) '-(142-147) convection,tot:d<uluiuj>/dxl'
      write(*,*) '-(148-153) convection,mean:d(<ul><ui><uj>)/dxl'
      write(*,*) "-(154-159) convection,fluct:<ul>d<ui'uj'>/dxl"
      write(*,*) '-(160-165) molecular diffusion,tot:1/Re<d2uiuj/dxl2>'
      write(*,*) '-(166-171) molecular diffusion,mean:'
      write(*,*) '-(172-177) molecular diffusion,fluct:'
      write(*,*) '-(178-183) turbulent diffusion mean: '
      write(*,*) "-(184-189) turbulentdiffusionfluct:-d<ui'uj'ul'>/dxl"
      write(*,*) "-(190-195) production -<ui'ul'>d<uj>/dxl"
      write(*,*) '-(196-201) v-p_gra,tot:-<uidp/dxj-ujdp/dxi>'
      write(*,*) '-(202-207) v-p_gra,mean:-<ui>d<p>/dxj-<uj>d<p>/dxi'
      write(*,*) "-(208-213) v-p_gra,fluct: -<ui'dp'/dxj>-<uj'dp'/dxi>"
      write(*,*) '-(214-218) (v-p)_gra,tot: -<dpui/dxj>-<dpuj/dxi>'
      write(*,*) '-(219-223) (v-p)_gra,mean:'
      write(*,*) "-(224-228) (v-p)_gra,fluct:-<dp'ui'/dxj>-<dp'uj'/dxi>"
      write(*,*) '-(229-234) p-strain rate,tot:2<p S_ij>'
      write(*,*) '-(235-240) p-strain rate,mean:2<p>< S_ij>'
      write(*,*) "-(241-246) p-strain rate,fluct:2<p'S_ij'> "
      write(*,*) '-(247-252) dissipation eps,tot:-2/Re<dui/dxl duj/dxl>'
      write(*,*) '-(253-258) dissipation eps,mean'
      write(*,*) '-(259-264) dissipation eps,fluct'
      write(*,*) '-(265-267) SGS diffusion, tot, mean, fluct'
      write(*,*) '-(268) total shear stress'
      write(*,*) '-(269-272) mean/fluct Joule diss/MHD diffusion'
      write(*,*) '-(273-274) total budget mean/fluct'
      if (scalar.gt.0) then
         write(*,*) '0 for the statistics of  the scalar'
      end if
      read(*,*) ivar

c     AM
c      if (ivar.eq.-71) then
c         write(*,*) 'Choose if you want the first maximum or the global
c     & maximum. Choose 0 for the first maximum and 1 for the
c     & global maximum'
c         read(*,*)firstmax
c      end if
      if (ivar.eq.0) then
         if (scalar.eq.0) stop
      end if
      if (ivar.gt.nxys) stop
      if (ivar.lt.-nxys2) stop
      if (ivar.ne.0) goto 2000

      write(*,*) '1     theta'
      write(*,*) '2     theta_rms'
      write(*,*) '3-5   utheta, vtheta, wtheta'
      write(*,*) '6-8   utheta^2, vtheta^2, wtheta^2'
      write(*,*) '9-11  dtheta/dx_i dtheta/dx_i'
      write(*,*) '12    p theta'
      write(*,*) '13-15 p thetax, p thetay, p thetaz'
      write(*,*) '16-18 u^2 theta, v^2 theta, w^2 theta'
      write(*,*) '19-21 uvtheta, uwtheta, vwtheta'
      write(*,*) '22-24 u thetax, u thetay, u thetaz'
      write(*,*) '25-27 v thetax, v thetay, v thetaz'
      write(*,*) '28-30 w thetax, w thetay, w thetaz'
      write(*,*) '31-33 dui/dxj dtheta/dxj'
      write(*,*) '34-35 theta^3,theta^4'
      write(*,*) '36    sig_i dtheta/dx_i'
      write(*,*) '37-38 sig_i dtheta/dx_i (forward/backscatter)'
      write(*,*) '39-41 sig_i theta'
      write(*,*) '42-44 sig_i'
      write(*,*)
      write(*,*) '-(1-2)   thetax,thetay'
      write(*,*) '-(3-5)   <utheta>x,<utheta>y,<vtheta>x'
      write(*,*) '-(6-8)   <vtheta>y,<wtheta>x,<wtheta>y'
      write(*,*) '-(9-11)  theta ux,theta uy,theta uz'
      write(*,*) '-(12-14) theta vx,theta vy,theta vz'
      write(*,*) '-(15-17) theta wx,theta wy,theta wz'
      write(*,*) '-(18-19) skewness and flatness for the scalar'
      write(*,*) '-20      turbulent Pr'
      write(*,*) '-21      k-theta'
      write(*,*) '-(22-23) (k-theta)x,(k-theta)y'
      write(*,*) '-(24-25) (k-theta)xx,(k-theta)yy'
      write(*,*) '-26      th-con,tot: 1/2d<ththuj>/dxj'
      write(*,*) '-27      th-con,mean: 1/2d(<th><th><uj>)/dxj'
      write(*,*) "-28      th-con,fluct: 1/2<uj>d<th'th'>)/dxj"
      write(*,*) '-29      th-mol_dif,tot:1/(2Pe)*<d2th^2/dxjdxj>'
      write(*,*) '-30      th-mol_dif,mean: 1/(2Pe)*d2<th>^2/dxjdxj'
      write(*,*) "-31      th-mol_dif,fluct: 1/(2Pe)*<d2th'^2/dxjdxj>"
      write(*,*) "-32      th-tur_dif,mean: -d(<uj'th'><th>)/dxj"
      write(*,*) "-33      th-tur_dif,fluct: -1/2*d(<th'th'uj'>)/dxj"
      write(*,*) "-34      production -<ui'th'><dth/dxi>"
      write(*,*) '-35      eps_theta,tot:-1/Pe<dth/dxi dth/dxi>'
      write(*,*) '-36      eps_theta,mean:-1/Pe<dth/dxi><dth/dxi>'
      write(*,*) "-37      eps_theta,fluct:-1/Pe<dth'/dxi dth'/dxi>"
      write(*,*) '-(38-40) th-convection,tot: d<thuiuj>/dxj'
      write(*,*) '-(41-43) th-convection,mean: d(<th><ui><uj>)/dxj'
      write(*,*) "-(44-46) th-convection,fluct: <uj>d<ui'th'>)/dxj"
      write(*,*) '-(47-49) th-v mol_dif,tot:'
      write(*,*) '-(50-52) th-v mol_dif,mean:'
      write(*,*) '-(53-55) th-v mol_dif,fluct:'
      write(*,*) '-(56-58) th-v tur_dif,mean:'
      write(*,*) "-(59-61) th-v tur_dif,fluct: -d<ui'uj'th'>/dxj"
      write(*,*) '-(62-64) production'
      write(*,*) '-(65-66) (th-p)_gra,tot: -d<pth>/dxi'
      write(*,*) '-(67-68) (th-p)_gra,mean: -d<p><th>/dxi'
      write(*,*) "-(69-70) (th-p)_gra,fluct:-d<p'th'>/dxi"
      write(*,*) '-(71-73) p-th_gra,tot: <p dth/dxi>'
      write(*,*) '-(74-76) p-th_gra,mean: <p><dth/dxi>'
      write(*,*) '-(77-79) th-p_gra,tot: -<thdp/dxi>'
      write(*,*) '-(80-82) th-p_gra,mean: -<th>d<p>/dxi'
      write(*,*) "-(83-85) th-p_gra,fluct: -<th'dp'/dxi>"
      write(*,*) '-(86-88) th-v dis,tot:-(1/Re+1/Pe)<dth/dxj dui/dxj>'
      write(*,*) '-(89-91) th-v dis,mean:'
      write(*,*) '-(92-94) th-v dis,fluct:'
      write(*,*) '-95      eps_sgs_theta,tot: <sig_i dtheta/dx_i>'
      write(*,*) '-96      eps_sgs_theta,mean: <sig_i><dtheta/dx_i>'
      write(*,*) "-97      eps_sgs_theta,fluct: <sig_i' dtheta/dx_i'>"
      write(*,*) '-(98-100)sgs_theta_diffusion: tot, mean, fluct'
      read(*,*) ivarth

      if (ivarth.eq.0) stop
      if (ivarth.gt.nxysth) stop
      if (ivarth.lt.-nxysth2) stop
c
c     Default options
c
 2000 iopt=0
      mbox=1
      iplot=1
      ifilt=0
      iscale=0
      ihard=0
      lpgmr=.false.
      ltab=.false.
      lmat=.false.
      write(headb,9020) xl,zl,nx,ny,nz,namxys
 9020 format(' xl = ',F6.2,' zl= ',F6.2,' ',I3,'x',I3,'x',I3,
     &     ' file ',A20,'$')
c
c     Shifting to follow disturbance in poiseuille flow
c     xr is the coordinate of the center of the plot
c     xs is the coordinate in the center of the data (always 0.)
c
      xr=xl/2.
      dx=xl/nx
      mxr=int((xr-xs)/dx+.5)
      if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
      xr=xs+mxr*dx
c
c     Set optional parameters
c
      write(*,*) 'change default options ? (y/n)'
      read(*,'(a)') ans
      if (ans.ne.'N'.and.ans.ne.'n') then
        iopt=1
        write(*,*) 'type of output'
        write(*,*) '0 screen only'
        write(*,*) '1 postscript and screen'
        write(*,*) '2 postscript only'
        write(*,*) '3 tektronix file'
        write(*,*) '4 screen and pgmr file'
        write(*,*) '7 tek file and pgmr file'
        write(*,*) '8 screen and matlab file'
        write(*,*) '10 screen and table'
        read(*,*) ihard
        if (ihard.eq.8) then
          write(*,*) 'give file to write to'
          read(*,'(a)') xmat
          lmat=.true.
          ihard=ihard-8
        end if
        if (ihard.ge.10) then
          ltab=.true.
          ihard=ihard-10
        end if
        if (ihard.ge.4) then
          lpgmr=.true.
          ihard=ihard-4
        end if
c
c     Filter
c
        write(*,*) 'choose filter'
        write(*,*) '0 no filter'
        write(*,*) '1 gaussian lowpass filter'
        read(*,*) ifilt
        if (ifilt.ge.1) then
          write(*,*) 'give filter lengthscale l'
          write(*,*) 'the filter is exp(-(k*l/2pi)**2)'
          read(*,*) sfilt
        end if
c
c     Aspect ratio
c
        write(*,*) 'give type of plot'
        write(*,*) '1 xy contour plot'
        write(*,*) '2 f(x) at multiple y'
        write(*,*) '3 f(y) at multiple x'
        read(*,*) iplot
        if (iplot.eq.2) then
 1010     write(*,*) 'give number of y-positions, max ',mgr
          read(*,*) npy
          if (npy.gt.mgr) goto 1010
          do 2010 i=1,npy
            write(*,9030) i
 9030       format('give y-value for curve',i4)
            read(*,*) yp(i)
            if (fltype.eq.1.or.fltype.eq.2.or.
     &           fltype.eq.4.or.fltype.eq.5) then
              y1(i)=int(acos(yp(i))/pi*real(nyp-1)+1.5)
            else
              y1(i)=int(acos(yp(i)*dstar-1.)/pi*real(nyp-1)+1.5)
            end if
            yp(i)=eta(y1(i))
            y1(i)=nyp+1-y1(i)
            write(*,9040) yp(i),y1(i)
 9040       format(' nearest plane y=',F6.3,' plane number',I4)
 2010     continue
        end if
        if (iplot.eq.3) then
 1020     write(*,*) 'give number of x-positions, max ',mgr
          read(*,*) npx
          if (npx.gt.mgr) goto 1020
          do 2020 i=1,npx
            write(*,9050) i
 9050       format('give x-value for curve',i4)
            read(*,*) xp(i)
            if (xp(i).lt.0) then
              x1(i)=nx-mod(nx-1+int(-xp(i)/xl*real(nx)+.5),nx)
            else
              x1(i)=mod(int(xp(i)/xl*real(nx)+.5),nx)+1
            end if
            xp(i)=xl/real(nx)*real(x1(i)-1)
            write(*,9060) xp(i),x1(i)
 9060       format(' nearest plane x=',F8.3,' plane number',I4)
 2020     continue
        end if
        write(*,*) 'give scaling'
        write(*,*) '0 scale with inflow freestream vel and displ thick'
        write(*,*) '1 scale with local freestream vel and displ thick'
        write(*,*) '2 scale with local wall units'
        if (iplot.eq.1.or.iplot.eq.3) then
          write(*,*) '3 same as 1, and y scaled with local displ thick'
          write(*,*) '4 same as 2, and y scaled with local wall units'
          write(*,*) '5 same as 4, and divided by y+'
        end if
        read(*,*) iscale
        write(*,*) 'choose aspect ratio'
        write(*,*) '0 equal scale for vertical and horizontal'
        write(*,*) '1 largest possible picture (aspect ratio 1.6)'
        read(*,*) mbox
c
c     Boxshift
c
        if (iplot.le.2) then
          write(*,*) 'the center of the plot is now at x=',xr
          write(*,*) 'give offset (0 for no change)'
          read(*,*) xrd
          xr=xr+xrd
          mxr=int((xr-xs)/dx+.5)
          if (xr.lt.xs) mxr=-int((xs-xr)/dx+.5)
          xr=xs+mxr*dx
        end if
      else
c
c     Default options for special plots
c
         if (ivar.eq.-71.or.ivar.eq.-72) then
            iopt=0
            mbox=1
            iplot=2
            npy=1
            yp(1) = eta(15)
            y1(1)=nyp+1.-15.
            ifilt=0
            iscale=0
            ihard=0
            lpgmr=.false.
            ltab=.false.
            lmat=.true.
            xmat='urms'
         end if
      end if
      if (iplot.eq.1) then
        write(*,*) 'number of contours ? (0 to select manually)'
        read(*,*) ncont
        ncont=-ncont
        if (ncont.eq.0) then
          write(*,*)  'give start,end,spacing for contour levels'
          read(*,*) blev,elev,dlev
          ncont=int((elev-blev)/dlev+1.5)
          elev=blev+dlev*real(ncont-1)
          write(*,*) 'endlevel ',elev
        end if
      end if
c
c     Set computed parameters
c
      mx=nx+1
      my=nyp
      if (ivar.ge.19.and.ivar.le.22) my=nyp-1
      if (ivar.eq.1) var='      u'
      if (ivar.eq.2) var='      v'
      if (ivar.eq.3) var='      w'
      if (ivar.eq.4) var='   urms'
      if (ivar.eq.5) var='   vrms'
      if (ivar.eq.6) var='   wrms'
      if (ivar.eq.7) var='    omx'
      if (ivar.eq.8) var='    omy'
      if (ivar.eq.9) var='    omz'
      if (ivar.eq.10) var=' omxrms'
      if (ivar.eq.11) var=' omyrms'
      if (ivar.eq.12) var=' omzrms'
      if (ivar.eq.13) var='     uv'
      if (ivar.eq.14) var='     uw'
      if (ivar.eq.15) var='     vw'
      if (ivar.eq.16) var='u x*x+1'
      if (ivar.eq.17) var='v x*x+1'
      if (ivar.eq.18) var='w x*x+1'
      if (ivar.eq.19) var='u y*y+1'
      if (ivar.eq.20) var='v y*y+1'
      if (ivar.eq.21) var='w y*y+1'
      if (ivar.eq.22) var='u z*z+1'
      if (ivar.eq.23) var='v z*z+1'
      if (ivar.eq.24) var='w z*z+1'
      if (ivar.eq.25) var='reeps11'
      if (ivar.eq.26) var='reeps22'
      if (ivar.eq.27) var='reeps33'
      if (ivar.eq.28) var='reeps12'
      if (ivar.eq.29) var='reeps13'
      if (ivar.eq.30) var='reeps23'
      if (ivar.eq.31) var='pressure'
      if (ivar.eq.32) var='prms'
      if (ivar.eq.33) var='pu'
      if (ivar.eq.34) var='pv'
      if (ivar.eq.35) var='pw'
      if (ivar.eq.36) var='pux'
      if (ivar.eq.37) var='pvy'
      if (ivar.eq.38) var='pwz'
      if (ivar.eq.39) var='puy'
      if (ivar.eq.40) var='pvz'
      if (ivar.eq.41) var='pwx'
      if (ivar.eq.42) var='puz'
      if (ivar.eq.43) var='min u'
      if (ivar.eq.44) var='max u'
      if (ivar.eq.45) var='SijSij'
      if (ivar.eq.46) var='S11'
      if (ivar.eq.47) var='S12'
      if (ivar.eq.48) var='S13'
      if (ivar.eq.49) var='S22'
      if (ivar.eq.50) var='S23'
      if (ivar.eq.51) var='S33'
      if (ivar.eq.52) var='C=CS^2'
      if (ivar.eq.53) var='tijSij'
      if (ivar.eq.54) var='tau11'
      if (ivar.eq.55) var='tau12'
      if (ivar.eq.56) var='tau13'
      if (ivar.eq.57) var='tau22'
      if (ivar.eq.58) var='tau23'
      if (ivar.eq.59) var='tau33'
      if (ivar.eq.60) var='nu_t'
      if (ivar.eq.61) var='tijSij+'
      if (ivar.eq.62) var='tijSij-'
      if (ivar.eq.63) var='u^3'
      if (ivar.eq.64) var='v^3'
      if (ivar.eq.65) var='w^3'
      if (ivar.eq.66) var='v u^2'
      if (ivar.eq.67) var='w u^2'
      if (ivar.eq.68) var='u v^2'
      if (ivar.eq.69) var='w v^2'
      if (ivar.eq.70) var='u w^2'
      if (ivar.eq.71) var='v w^2'
      if (ivar.eq.72) var='uvw'
      if (ivar.eq.73) var='pvx'
      if (ivar.eq.74) var='pwy'
      if (ivar.eq.75) var='u^4'
      if (ivar.eq.76) var='v^4'
      if (ivar.eq.77) var='v^4'
      if (ivar.eq.78) var='p^3'
      if (ivar.eq.79) var='p^4'
      if (ivar.eq.80) var='uj tau1j'
      if (ivar.eq.81) var='uj tau2j'
      if (ivar.eq.82) var='uj tau3j'
      if (ivar.eq.83) var='phi    '
      if (ivar.eq.84) var='phi_rms'
      if (ivar.eq.85) var='j1     '
      if (ivar.eq.86) var='j2     '
      if (ivar.eq.87) var='j3     '
      if (ivar.eq.88) var='j1_rms '
      if (ivar.eq.89) var='j2_rms '
      if (ivar.eq.90) var='j3_rms '
      if (ivar.eq.91) var='<j1j2> '
      if (ivar.eq.92) var='<j1j3> '
      if (ivar.eq.93) var='<j2j3> '
      if (ivar.eq.94) var='phi j1 '
      if (ivar.eq.95) var='phi j2 '
      if (ivar.eq.96) var='phi j3 '

      if (ivar.eq.-1) var='      k'
      if (ivar.eq.-2) var=' q2ro11'
      if (ivar.eq.-3) var=' q2ro22'
      if (ivar.eq.-4) var=' q2ro33'
      if (ivar.eq.-5) var='lpurms2'
      if (ivar.eq.-6) var='lpvrms2'
      if (ivar.eq.-7) var='lpwrms2'
      if (ivar.eq.-8) var='lambdac'
      if (ivar.eq.-9) var='re_lamc'
      if (ivar.eq.-10) var=' eps11c'
      if (ivar.eq.-11) var=' eps22c'
      if (ivar.eq.-12) var=' eps33c'
      if (ivar.eq.-13) var='epsiloc'
      if (ivar.eq.-14) var='ux     '
      if (ivar.eq.-15) var='vx     '
      if (ivar.eq.-16) var='wx     '
      if (ivar.eq.-17) var='uy     '
      if (ivar.eq.-18) var='vy     '
      if (ivar.eq.-19) var='wy     '
      if (ivar.eq.-20) var='  eps11'
      if (ivar.eq.-21) var='  eps22'
      if (ivar.eq.-22) var='  eps33'
      if (ivar.eq.-23) var='  eps12'
      if (ivar.eq.-24) var='  eps13'
      if (ivar.eq.-25) var='  eps23'
      if (ivar.eq.-26) var='epsilon'
      if (ivar.eq.-27) var=' -uv/2k'
      if (ivar.eq.-28) var=' lambda'
      if (ivar.eq.-29) var='kx     '
      if (ivar.eq.-30) var='epsx   '
      if (ivar.eq.-31) var='ky     '
      if (ivar.eq.-32) var='epsy   '
      if (ivar.eq.-33) var='(pu)x'
      if (ivar.eq.-34) var='(pv)x'
      if (ivar.eq.-35) var='(pw)x'
      if (ivar.eq.-36) var='(pu)y'
      if (ivar.eq.-37) var='(pv)y'
      if (ivar.eq.-38) var='(pw)y'
      if (ivar.eq.-39) var='upx'
      if (ivar.eq.-40) var='vpx'
      if (ivar.eq.-41) var='upy'
      if (ivar.eq.-42) var='vpy'
      if (ivar.eq.-47) var='urms^2xx'
      if (ivar.eq.-48) var='urms^2yy'
      if (ivar.eq.-49) var='vrms^2xx'
      if (ivar.eq.-50) var='vrms^2yy'
      if (ivar.eq.-51) var='(uv)xx  '
      if (ivar.eq.-52) var='(uv)yy  '
      if (ivar.eq.-53) var='urms^2x '
      if (ivar.eq.-54) var='vrms^2y '
      if (ivar.eq.-55) var='urms^2x '
      if (ivar.eq.-56) var='vrms^2y '
      if (ivar.eq.-57) var='(uv)x   '
      if (ivar.eq.-58) var='(uv)y   '
      if (ivar.eq.-59) var='kyy     '
      if (ivar.eq.-60) var='kxx     '
      if (ivar.eq.-61) var='uxx     '
      if (ivar.eq.-62) var='vxx     '
      if (ivar.eq.-63) var='wxx     '
      if (ivar.eq.-64) var='uyy     '
      if (ivar.eq.-65) var='vyy     '
      if (ivar.eq.-66) var='wyy     '
      if (ivar.eq.-67) var='Lk      '
      if (ivar.eq.-68) var='Lu      '
      if (ivar.eq.-69) var='Lu / Lk '
      if (ivar.eq.-70) var='Lu/Lk n.'
      if (ivar.eq.-71) var='m(urms) '
      if (ivar.eq.-72) var='ym(urms)'
      if (ivar.eq.-73) var='Tu(lok) '
      if (ivar.eq.-74) var='Tu(free)'
      if (ivar.eq.-75) var='Blas    '
      if (ivar.eq.-76) var='dev     '
      if (ivar.eq.-77) var='dev max '
      if (ivar.eq.-78) var='dev min '
      if (ivar.eq.-82) var='evtot   '
      if (ivar.eq.-83) var='evmean  '
      if (ivar.eq.-84) var='evfluct '
      if (ivar.eq.-85) var='C_S     '
      if (ivar.eq.-86) var='estot   '
      if (ivar.eq.-87) var='esmean  '
      if (ivar.eq.-88) var='esfluct '
      if (ivar.eq.-89) var='P       '
      if (ivar.eq.-90) var='P/eps   '
      if (ivar.eq.-91) var='u1u1    '
      if (ivar.eq.-92) var='u2u2    '
      if (ivar.eq.-93) var='u3u3    '
      if (ivar.eq.-94) var='u1u2    '
      if (ivar.eq.-99) var='wrms^2yy'
      if (ivar.eq.-100) var='wrms^2x '
      if (ivar.eq.-101) var='wrms^2y '
      if (ivar.eq.-102) var='(uw)xx  '
      if (ivar.eq.-103) var='(uw)yy  '
      if (ivar.eq.-104) var='(vw)xx  '
      if (ivar.eq.-105) var='(vw)yy  '
      if (ivar.eq.-106) var='(uw)x   '
      if (ivar.eq.-107) var='(uw)y   '
      if (ivar.eq.-108) var='(vw)x   '
      if (ivar.eq.-109) var='(vw)y   '
      if (ivar.eq.-110) var='wpx     '
      if (ivar.eq.-111) var='wpy     '
      if (ivar.eq.-112) var='upz     '
      if (ivar.eq.-113) var='vpz     '
      if (ivar.eq.-114) var='wpz     '
      if (ivar.eq.-115) var='px      '
      if (ivar.eq.-116) var='py      '
      if (ivar.eq.-117) var='s(u)    '
      if (ivar.eq.-118) var='s(v)    '
      if (ivar.eq.-119) var='s(w)    '
      if (ivar.eq.-120) var='f(u)    '
      if (ivar.eq.-121) var='f(v)    '
      if (ivar.eq.-122) var='f(w)    '
      if (ivar.eq.-123) var='s(p)    '
      if (ivar.eq.-124) var='f(p)    '

      if (ivar.eq.-265) var='sdiftot'
      if (ivar.eq.-266) var='sdifmean'
      if (ivar.eq.-267) var='sdiffluc'
      if (ivar.eq.-268) var='totshear'
      if (ivar.eq.-269) var='mhd_dism'
      if (ivar.eq.-270) var='mhd_disf'
      if (ivar.eq.-271) var='mhd_difm'
      if (ivar.eq.-272) var='mhd_diff'
      if (ivar.eq.-273) var='tot budm'
      if (ivar.eq.-274) var='tot budf'

      if (ivarth.eq.1) var='  theta'
      if (ivarth.eq.2) var=' th_rms'
      if (ivarth.eq.3) var='u theta'
      if (ivarth.eq.4) var='v theta'
      if (ivarth.eq.5) var='w theta'
      if (ivarth.eq.6) var='utheta^2'
      if (ivarth.eq.7) var='vtheta^2'
      if (ivarth.eq.8) var='wtheta^2'
      if (ivarth.eq.9) var='th_x^2 '
      if (ivarth.eq.10) var='th_y^2'
      if (ivarth.eq.11) var='th_z^2'
      if (ivarth.eq.12) var='p theta'
      if (ivarth.eq.13) var='pthetax'
      if (ivarth.eq.14) var='pthetay'
      if (ivarth.eq.15) var='pthetaz'
      if (ivarth.eq.16) var='u^2 th '
      if (ivarth.eq.17) var='v^2 th '
      if (ivarth.eq.18) var='w^2 th '
      if (ivarth.eq.19) var='uvtheta'
      if (ivarth.eq.20) var='uwtheta'
      if (ivarth.eq.21) var='vwtheta'
      if (ivarth.eq.22) var='uthetax'
      if (ivarth.eq.23) var='uthetay'
      if (ivarth.eq.24) var='uthetaz'
      if (ivarth.eq.25) var='vthetax'
      if (ivarth.eq.26) var='vthetay'
      if (ivarth.eq.27) var='vthetaz'
      if (ivarth.eq.28) var='wthetax'
      if (ivarth.eq.29) var='wthetay'
      if (ivarth.eq.30) var='wthetaz'
      if (ivarth.eq.31) var='u_jth_j'
      if (ivarth.eq.32) var='v_jth_j'
      if (ivarth.eq.33) var='w_jth_j'
      if (ivarth.eq.34) var='theta^3'
      if (ivarth.eq.35) var='theta^4'
      if (ivarth.eq.36) var=' sidthi'
      if (ivarth.eq.37) var='sidthi+'
      if (ivarth.eq.38) var='sidthi-'
      if (ivarth.eq.39) var='sig1 th'
      if (ivarth.eq.40) var='sig2 th'
      if (ivarth.eq.41) var='sig3 th'
      if (ivarth.eq.42) var='sig_1'
      if (ivarth.eq.43) var='sig_2'
      if (ivarth.eq.44) var='sig_3'

      if (ivarth.eq.-1) var='thetax '
      if (ivarth.eq.-2) var='thetay '
      if (ivarth.eq.-3) var='(uth)x '
      if (ivarth.eq.-4) var='(uth)y '
      if (ivarth.eq.-5) var='(vth)x '
      if (ivarth.eq.-6) var='(vth)y '
      if (ivarth.eq.-7) var='(wth)x '
      if (ivarth.eq.-8) var='(wth)y '
      if (ivarth.eq.-9) var='th ux  '
      if (ivarth.eq.-10) var='th uy  '
      if (ivarth.eq.-11) var='th uz  '
      if (ivarth.eq.-12) var='th vx  '
      if (ivarth.eq.-13) var='th vy  '
      if (ivarth.eq.-14) var='th vz  '
      if (ivarth.eq.-15) var='th wx  '
      if (ivarth.eq.-16) var='th wy  '
      if (ivarth.eq.-17) var='th wz  '
      if (ivarth.eq.-18) var='s(th)  '
      if (ivarth.eq.-19) var='f(th)  '
      if (ivarth.eq.-20) var='Pr_t   '
      if (ivarth.eq.-21) var='k_theta'
      if (ivarth.eq.-22) var='(k_th)x'
      if (ivarth.eq.-23) var='(k_th)y'
      if (ivarth.eq.-24) var='k_th xx'
      if (ivarth.eq.-25) var='k_th yy'
      if (ivarth.eq.-95) var='esttot'
      if (ivarth.eq.-96) var='estmean'
      if (ivarth.eq.-97) var='estfluc'
      if (ivarth.eq.-98) var='sdiftot'
      if (ivarth.eq.-99) var='sdifmean'
      if (ivarth.eq.-100) var='sdiffluct'
      scl=' '
      if (iscale.eq.2.or.iscale.eq.4) scl='+'
      if (iscale.eq.1.or.iscale.eq.3) scl=' outer scale'
      if (iscale.eq.5) scl='+/y+'
c
c     Fetch the selected variable and make a grid
c
      call mgrid(plane,pxys,pxys2,pxysth,pxysth2,nxys2,nxysth2,
     &     ivar,ivarth,mx,my,grid,xr,xl,mxr,eta)
c
c     Scale the independent and dependent variable
c
      if (iscale.gt.0) call rscale(plane,grid,mxr,mx,my,ivar,ivarth,
     & iscale,uinf,utau,delta,nux,pdim,pdim2,thdim,thdim2,nxys,nxys2,
     &     nxysth,nxysth2,re,thdelta,tinf,twall)
c
c     Low-pass filter the dependent variable
c
      if (ifilt.ge.1) call filter(plane,mx,my,ifilt,sfilt,xl)
c
      write(heada,9510) var,scl,t,re
 9510   format(1x,A7,A12,' t= ',F6.1,' re= ',F7.1,'$')
      if (iplot.eq.1) then
         call cont5(mx,my,plane,grid,iopt,
     &        ncont,blev,elev,'x$','y$',heada,headb,mbox,ihard,
     &       1,5,6,lpgmr)
         if (lmat) then
            write(*,*) 'write to ',xmat
            open(unit=24,file=xmat,form='unformatted')
            write(24) mx,my
            write(24) ((grid(x+(y-1)*mx+(1-1)*mx*my,1,1),
     &           x=1,mx),y=1,my)
            write(24) ((grid(x+(y-1)*mx+(2-1)*mx*my,1,1),
     &           x=1,mx),y=1,my)
            write(24) ((plane(x+(y-1)*mx,1),x=1,mx),y=1,my)
            close(24)
         end if

      end if
      if (iplot.eq.2) then
        if (npy.eq.1) then
          write(headc,9520) (yp(i),i=1,1)
 9520     format(' y= ',F8.3,'$')
        end if
        if (npy.eq.2) then
          write(headc,9521) (yp(i),i=1,2)
 9521     format(' solid: y= ',F8.3,' dashed: y= ',F8.3,'$')
        end if
        if (npy.eq.3) then
          write(headc,9522) (yp(i),i=1,3)
 9522     format(' solid: y= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,'$')
        end if
        if (npy.eq.4) then
          write(headc,9523) (yp(i),i=1,4)
 9523     format(' solid: y= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,' chdh: ',f8.3,'$')
        end if
        if (npy.eq.5) then
          write(headc,9524) (yp(i),i=1,5)
 9524     format(' solid: y= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,' chdh: ',f8.3,' chdt: ', f8.3,'$')
        end if
        do 2100 i=1,npy
          npoint(i)=mx
 2100   continue
        call mxpl(xgraph,xgrid,plane,grid,y1,npy,mx,my,mgr)
        if (ltab) then
          write(49,*) npy
          do 2102 i=1,npy
          write(49,*) yp(i),mx
          do 2102 x=1,mx
            write(49,*)x,xgrid(x,i),xgraph(x,i)
 2102     continue
        end if
        if (lmat) then
          open(unit=24,file=xmat)
c          do 2103 i=1,npy
          do 2103 y=1,mx
             if (npy.eq.5) then
            write(24,'(6e18.9)')xgrid(y,1),(xgraph(y,i),i=1,npy)
            else
            write(24,*)xgrid(y,1),(xgraph(y,i),i=1,npy)
            end if
 2103     continue
        end if
        call rita1a(xgrid,xgraph,0.,0.,0.,0.,npoint,npy,mx,
     &       'x',var,heada,headb,headc,mbox,1,ihard)
      end if
      if (iplot.eq.3) then
        if (npx.eq.1) then
          write(headc,9530) (xp(i),i=1,1)
 9530     format(' x= ',F8.3,'$')
        end if
        if (npx.eq.2) then
          write(headc,9531) (xp(i),i=1,2)
 9531     format(' solid: x= ',F8.3,' dashed: x= ',F8.3,'$')
        end if
        if (npx.eq.3) then
          write(headc,9532) (xp(i),i=1,3)
 9532     format(' solid: x= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,'$')
        end if
        if (npx.eq.4) then
          write(headc,9533) (xp(i),i=1,4)
 9533     format(' solid: x= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,' chdh: ',f8.3,'$')
        end if
        if (npx.eq.5) then
          write(headc,9534) (xp(i),i=1,5)
 9534     format(' solid: x= ',F8.3,' dash: ',F8.3,
     &    ' dot: ',f8.3,' chdh: ',f8.3,' chdt: ',f8.3,'$')
        end if
        do 2110 i=1,npx
          npoint(i)=my
 2110   continue
        call mypl(ygraph,ygrid,plane,grid,x1,npx,mx,my,mgr)
        xscale='y'
        if (iscale.eq.3) xscale='y/dstar'
        if (iscale.eq.4.or.iscale.eq.5) xscale='y+'
        if (ltab) then
          write(49,*) npx
          do 2105 i=1,npx
            write(49,*) xp(i),my
          do 2105 y=1,my
            write(49,*)y,ygrid(y,i),ygraph(y,i)
 2105     continue
        end if
        if (lmat) then
c
c     Write Matlab etc. output
c
          open(unit=24,file=xmat)
          do y=1,my
             write(24,'(SP100ES26.16E3)') ygrid(y,1),
     &            (ygraph(y,i),i=1,npx)
          end do
          close(24)
        end if
c
c     Do the Tektronix plot
c
        call rita1a(ygrid,ygraph,0.,0.,0.,0.,npoint,npx,my,
     &       xscale,var,heada,headb,headc,mbox,1,ihard)

      end if
c
c     Next plot
c
      goto 1000

      end program pxyst
