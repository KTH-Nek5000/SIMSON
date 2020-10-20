c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine fring(om2r,om2i,u2r,u2i,tc,xsc,zsc,xl,zl,yb,
     &     fstart,fend,bu1,bu2,
     &     osmod,osnumb,osur,osui,osvr,osvi,oswr,oswi,
     &     evr,evi,eur,eui,ewr,ewi,afw,alfaw,betaw,ampob,amp2d,
     &     ev2r,ev2i,eu2r,eu2i,afw2d,alf2d,
     &     xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &     th2r,th2i,
     &     ampst,streak,betast,omegast,
     &     ndxst,uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,
     &     tsmoo,tsmst,tsmend,iampst,phist,pert)

      implicit none

      include 'par.f'

      integer yb
      real tc,xsc,zsc,xl,zl
      real fstart,fend
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real afw,alfaw,betaw,ampob,amp2d
      real ev2r(nyp),ev2i(nyp),eu2r(nyp),eu2i(nyp),afw2d,alf2d
      real evr(nyp),evi(nyp),eur(nyp),eui(nyp),ewr(nyp),ewi(nyp)
      real u2r(nxp/2+1,mby,nzd,3),u2i(nxp/2+1,mby,nzd,3)
      real om2r(nxp/2+1,mby,nzd,3),om2i(nxp/2+1,mby,nzd,3)
      real th2r(nxp/2+1,mby,nzd,4*scalar),th2i(nxp/2+1,mby,nzd,4*scalar)

      real xc1(nxp/2),xc2(nxp/2),csxc(nxp/2+1,8)
      real zc
      integer x,y,z,y1,xstart,xend,xendc,xst,xen,i,imin
      real c1s,c1c,c2s,c2c,c1sw,c1cw,c2sw,c2cw,bz1,bz2
      real k1s,k1c,k2s,k2c
      real fring1(nxp/2),fring2(nxp/2),afwtc,afw2dt
      real fstc,fenc
      real r
      parameter(r=0.499999999999999999)
      logical osmod
      integer osnumb
      real osur(osnf,nyp), osvr(osnf,nyp), oswr(osnf,nyp)
      real osui(osnf,nyp), osvi(osnf,nyp), oswi(osnf,nyp)
c
c     Streak generation
c
      logical streak
      real iampst,ampst(2),tsmoo(4),tsmst(2),tsmend(2)
      real betast(2),omegast(2),phist
      real uust_r(nyp,180,2),vvst_r(nyp,180,2),wwst_r(nyp,180,2)
      real uust_i(nyp,180,2),vvst_i(nyp,180,2),wwst_i(nyp,180,2)
      integer ndxst,xshift
      real smooth,stamps(2),bbz1,bbz2,tzero

      real cphas1(osnf,nxp/2,nzd),cphas2(osnf,nxp/2,nzd)
      real sphas1(osnf,nxp/2,nzd),sphas2(osnf,nxp/2,nzd)
      real udis1,udis2,vdis1,vdis2,wdis1,wdis2
      real ff,zll

      integer j,ith
      logical pert

      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Functions
c
      real,external :: step

c
c     Note we don't need to generate the fringe function every time step
c     unless the moving coordinate is used. Thus the fringe function can be
c     generated just once for the non-moving one in other routine, preprbl,
c     for the non-moving coordinate. Thus do not give fstart/fend across the
c     discontinuity. Here xst, xen are starting or ending points of non-zero
c     fringe region.
c
c     Evaluating the forcing takes a considerable CPU time in large cases
c     with ampob.gt.0. Therefore we only want to compute it if fringe is
c     non zero. If we are calculating in a moving frame of reference we have
c     to map fstart and fend to positions where they are inside the box
c
      fstc=fstart+int((xsc-fstart)/xl+r)*xl
      fenc=fend+int((xsc-fend)/xl+r)*xl
      xstart=int(real(nxp/2)*(fstc-xsc+xl/2.)/xl)+1
      xend=int(real(nxp/2)*(fenc-xsc+xl/2.)/xl)+1

      if (xstart.gt.xend) then
         imin=0
         xendc=nxp/2
      else
         imin=1
         xendc=1
      end if

      do i=imin,1
         xst=xstart**i
         xen=max(xend,xendc**i)

         if (streak.and.(ndxst.ne.(xen-xst+1)*2)) then
            write(*,*)'different number of points in streak file'
            write(*,*)' fields from input file= ',ndxst
            write(*,*)' fields from code= ',(xen-xst+1)*2
            write(*,*)' xsc      = ',xsc
            write(*,*)' xst, xen = ',xst,xen
c            call stopnow(45740)
            write(*,*) 'correcting...'
            xshift = (xen-xst+1)*2 - ndxst
            do j=1,2
               do x=180,1+xshift,-1
                  do y=1,nyp
                     uust_r(y,x,j)=uust_r(y,x-xshift,j)
                     uust_i(y,x,j)=uust_i(y,x-xshift,j)
                     vvst_r(y,x,j)=vvst_r(y,x-xshift,j)
                     vvst_i(y,x,j)=vvst_i(y,x-xshift,j)
                     wwst_r(y,x,j)=wwst_r(y,x-xshift,j)
                     wwst_i(y,x,j)=wwst_i(y,x-xshift,j)
                  end do
               end do
               do x=xshift,1,-1
                  do y=1,nyp
                     uust_r(y,x,j)=0.
                     uust_i(y,x,j)=0.
                     vvst_r(y,x,j)=0.
                     vvst_i(y,x,j)=0.
                     wwst_r(y,x,j)=0.
                     wwst_i(y,x,j)=0.
                  end do
               end do
            end do
            ndxst = (xen-xst+1)*2
         end if

         if (ampob.ne.0..or.amp2d.ne.0.) then
c
c     Damping to blasius flow plus forced oscillation in the fringe region.
c
c     Add one 2d wave in the fringe region and
c     Add two oblique waves in the fringe region. The generation
c     of the angles would in a strait forward manner be described by
c     a1=alfaw*xc1(x)+betaw*zc-afw*tc
c     b1=alfaw*xc1(x)-betaw*zc-afw*tc
c     a2=alfaw*xc2(x)+betaw*zc-afw*tc
c     b2=alfaw*xc2(x)-betaw*zc-afw*tc
c     c1c=cos(a1)+cos(b1)
c     c1s=sin(a1)+sin(b1)
c     c2c=cos(a2)+cos(b2)
c     c2s=sin(a2)+sin(b2)
c     c1cw=cos(a1)-cos(b1)
c     c1sw=sin(a1)-sin(b1)
c     c2cw=cos(a2)-cos(b2)
c     c2sw=sin(a2)-sin(b2)
c     this is done in a more efficient way using realtions for
c     sums and differences of sines and cosines
c
            afwtc=afw*tc
            afw2dt=afw2d*tc

            do x=xst,xen
               csxc(x,1)=cos(alfaw*xc1(x)-afwtc)
               csxc(x,2)=sin(alfaw*xc1(x)-afwtc)
               csxc(x,3)=cos(alfaw*xc2(x)-afwtc)
               csxc(x,4)=sin(alfaw*xc2(x)-afwtc)
               csxc(x,5)=cos(alf2d*xc1(x)-afw2dt)
               csxc(x,6)=sin(alf2d*xc1(x)-afw2dt)
               csxc(x,7)=cos(alf2d*xc2(x)-afw2dt)
               csxc(x,8)=sin(alf2d*xc2(x)-afw2dt)
            end do
            do z=1,nzpc
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
               bz1=2*cos(betaw*zc)
               bz2=2*sin(betaw*zc)
               do y=1,min(mby,nyp-yb+1)
                  y1=y+yb-1
                  do x=xst,xen
                     c1c=csxc(x,1)
                     c1s=csxc(x,2)
                     c2c=csxc(x,3)
                     c2s=csxc(x,4)
                     k1c=csxc(x,5)
                     k1s=csxc(x,6)
                     k2c=csxc(x,7)
                     k2s=csxc(x,8)
                     c1cw=-c1s*bz2
                     c1sw=c1c*bz2
                     c2cw=-c2s*bz2
                     c2sw=c2c*bz2
                     c1c=c1c*bz1
                     c1s=c1s*bz1
                     c2c=c2c*bz1
                     c2s=c2s*bz1
                     om2r(x,y,z,1)=om2r(x,y,z,1)-fring1(x)*(
     &                    -(eur(y1)*c1c-eui(y1)*c1s)
     &                    -(eu2r(y1)*k1c-eu2i(y1)*k1s))
                     om2i(x,y,z,1)=om2i(x,y,z,1)-fring2(x)*(
     &                    -(eur(y1)*c2c-eui(y1)*c2s)
     &                    -(eu2r(y1)*k2c-eu2i(y1)*k2s))
                     om2r(x,y,z,2)=om2r(x,y,z,2)-fring1(x)*(
     &                    -(evr(y1)*c1c-evi(y1)*c1s)
     &                    -(ev2r(y1)*k1c-ev2i(y1)*k1s))
                     om2i(x,y,z,2)=om2i(x,y,z,2)-fring2(x)*(
     &                    -(evr(y1)*c2c-evi(y1)*c2s)
     &                    -(ev2r(y1)*k2c-ev2i(y1)*k2s))
                     om2r(x,y,z,3)=om2r(x,y,z,3)-fring1(x)*(
     &                    -(ewr(y1)*c1cw-ewi(y1)*c1sw))
                     om2i(x,y,z,3)=om2i(x,y,z,3)-fring2(x)*(
     &                    -(ewr(y1)*c2cw-ewi(y1)*c2sw))
                  end do
               end do
            end do
         end if
c
c     Forcing of the Orr-Sommerfeld modes
c
c     Known parameters:
c     osbeta        (in z )
c     osomega       (in tc)
c     osalr, osali  (in x )
c     osur, osui
c     osvr, osvi    (in y )
c     oswr, oswi
c     in general we have for a=ar+I*ai, b=br+I*bi
c     real(a*exp(I*b)) = exp(-bi)*(ar*cos(br)-ai*sin(br))
c
c     Here we have
c     u=real(U*exp(I*(alpha*x+beta*y-omega*t)))
c     with alpha,U complex and beta,omega,x,y,t real
c
c     So we would have something like (for every component)
c          udis1   = osamp*real(cmplx(osur(y),osui(y))*
c     &              exp(cmplx(0.,1.)*(
c     &                 cmplx(osalr,osali)*xc1(x)
c     &              +  osbeta*zc
c     &              -  osomega*tc )))
c          udis2   = osamp*real(cmplx(osur(y),osui(y))*
c     &              exp(cmplx(0.,1.)*(
c     &                 cmplx(osalr,osali)*xc2(x)
c     &              +  osbeta*zc
c     &              -  osomega*tc )))
c
c     In order to make it efficient, use the above decomposition
c     u=exp(-bi)*(ar*cos(br)-ai*sin(br)) with
c     ar = osur (y)
c     ai = osui (y)
c     br = osalr*x + osbeta*z - osomega*t
c     bi = osali*x
c
c     Computation is now performed in ffun.f
c
         if (osmod) then
            do z=1,nzpc
               do y=1,min(mby,nyp-yb+1)
                  y1=y+yb-1
                  do x=xst,xen
                     udis1=0.
                     udis2=0.
                     vdis1=0.
                     vdis2=0.
                     wdis1=0.
                     wdis2=0.
                     do j=1,osnumb
                        udis1=udis1+osur(j,y1)*cphas1(j,x,z)-
     &                       osui(j,y1)*sphas1(j,x,z)
                        udis2=udis2+osur(j,y1)*cphas2(j,x,z)-
     &                       osui(j,y1)*sphas2(j,x,z)

                        vdis1=vdis1+osvr(j,y1)*cphas1(j,x,z)-
     &                       osvi(j,y1)*sphas1(j,x,z)
                        vdis2=vdis2+osvr(j,y1)*cphas2(j,x,z)-
     &                       osvi(j,y1)*sphas2(j,x,z)

                        wdis1=wdis1+oswr(j,y1)*cphas1(j,x,z)-
     &                       oswi(j,y1)*sphas1(j,x,z)
                        wdis2=wdis2+oswr(j,y1)*cphas2(j,x,z)-
     &                       oswi(j,y1)*sphas2(j,x,z)
                     end do
                     om2r(x,y,z,1)=om2r(x,y,z,1)-fring1(x)*
     &                    (- udis1)
                     om2i(x,y,z,1)=om2i(x,y,z,1)-fring2(x)*
     &                    (- udis2)
                     om2r(x,y,z,2)=om2r(x,y,z,2)-fring1(x)*
     &                    (- vdis1)
                     om2i(x,y,z,2)=om2i(x,y,z,2)-fring2(x)*
     &                    (- vdis2)
                     om2r(x,y,z,3)=om2r(x,y,z,3)-fring1(x)*
     &                    (- wdis1)
                     om2i(x,y,z,3)=om2i(x,y,z,3)-fring2(x)*
     &                    (- wdis2)
                  end do
               end do
            end do
         end if
c
c     No wave forcing, damping to Blasius flow in the fringe region
c
         if (.not.pert) then
            do z=1,nzpc
               do y=1,min(mby,nyp-yb+1)
                  y1=y+yb-1
                  do x=xst,xen
                     om2r(x,y,z,1)=om2r(x,y,z,1)-
     &                    fring1(x)*(u2r(x,y,z,1)-bu1(x,y1,1))
                     om2i(x,y,z,1)=om2i(x,y,z,1)-
     &                    fring2(x)*(u2i(x,y,z,1)-bu2(x,y1,1))
                     om2r(x,y,z,2)=om2r(x,y,z,2)-
     &                    fring1(x)*(u2r(x,y,z,2)-bu1(x,y1,2))
                     om2i(x,y,z,2)=om2i(x,y,z,2)-
     &                    fring2(x)*(u2i(x,y,z,2)-bu2(x,y1,2))
                     om2r(x,y,z,3)=om2r(x,y,z,3)-
     &                    fring1(x)*(u2r(x,y,z,3)-bu1(x,y1,3))
                     om2i(x,y,z,3)=om2i(x,y,z,3)-
     &                    fring2(x)*(u2i(x,y,z,3)-bu2(x,y1,3))
                  end do
               end do
            end do
         else
c
c     Damp to zero for perturbation mode
c
            do z=1,nzpc
               do y=1,min(mby,nyp-yb+1)
                  y1=y+yb-1
                  do x=xst,xen
                     om2r(x,y,z,1)=om2r(x,y,z,1)-
     &                    fring1(x)*(u2r(x,y,z,1))
                     om2i(x,y,z,1)=om2i(x,y,z,1)-
     &                    fring2(x)*(u2i(x,y,z,1))
                     om2r(x,y,z,2)=om2r(x,y,z,2)-
     &                    fring1(x)*(u2r(x,y,z,2))
                     om2i(x,y,z,2)=om2i(x,y,z,2)-
     &                    fring2(x)*(u2i(x,y,z,2))
                     om2r(x,y,z,3)=om2r(x,y,z,3)-
     &                    fring1(x)*(u2r(x,y,z,3))
                     om2i(x,y,z,3)=om2i(x,y,z,3)-
     &                    fring2(x)*(u2i(x,y,z,3))
                  end do
               end do
            end do
         end if
c
c     Force two streaks in the fringe region
c     Set amplitude according to smoothing
c
         if (streak) then
            smooth=step( (tc-tsmst (1)) / tsmoo(1) )
     &           - step( (tc-tsmend(1)) / tsmoo(2) )
            stamps(1)=iampst+smooth*(ampst(1)-iampst)


c     To have a linear variation of the amplitude
c            stamps(1) = ampst(1)*(tc-tsmst(1))/tsmoo(1)
c            stamps(1) = min(max(stamps(1),0.),ampst(1))
c           write(*,*) stamps(1),tc,tsmst(1),tsmoo(1)

            smooth=step((tc-tsmst(2))/tsmoo(3))
     &           -step((tc-tsmend(2))/tsmoo(4))
            stamps(2)=iampst+smooth*(ampst(2)-iampst)

            tzero=(tsmend(1)+tsmst(2)+tsmoo(3))/2.

c     To completely kill streak 2
c            tzero = 0.
c            stamps(2) = 0.

            zll = 1.0/0.3
            do z=1,nzpc
c
c     Cosine/sine in z
c
               zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
c
c     Only half of the domain has streaks
c               ff = (1.-step(zc*zll))*step((zc+zl/2)*zll)
c     Otherwise
               ff=1.

               bz1=cos(betast(1)*zc-omegast(1)*(tc-tzero))
               bz2=sin(betast(1)*zc-omegast(1)*(tc-tzero))

               bbz1=cos(betast(2)*zc-omegast(2)*(tc-tzero)+phist*pi)
               bbz2=sin(betast(2)*zc-omegast(2)*(tc-tzero)+phist*pi)

               do x=xst,xen
                  do y=1,min(mby,nyp-yb+1)
                     y1=y+yb-1
                     om2r(x,y,z,1)=om2r(x,y,z,1)-fring1(x)*(
     &                    -stamps(1)*ff*
     &                    (uust_r(y1,(x-xst)*2+1,1)*bz1
     &                    -uust_i(y1,(x-xst)*2+1,1)*bz2)
     &                    -stamps(2)*ff*
     &                    (uust_r(y1,(x-xst)*2+1,2)*bbz1
     &                    -uust_i(y1,(x-xst)*2+1,2)*bbz2)
     &                    )

                     om2i(x,y,z,1)=om2i(x,y,z,1)-fring2(x)*(
     &                    -stamps(1)*ff*
     &                    (uust_r(y1,(x-xst+1)*2,1)*bz1
     &                    -uust_i(y1,(x-xst+1)*2,1)*bz2)-stamps(2)*ff*
     &                    (uust_r(y1,(x-xst+1)*2,2)*bbz1
     &                    -uust_i(y1,(x-xst+1)*2,2)*bbz2)
     &                    )

                     om2r(x,y,z,2)=om2r(x,y,z,2)-fring1(x)*(
     &                    -stamps(1)*ff*
     &                    (vvst_r(y1,(x-xst)*2+1,1)*bz1
     &                    -vvst_i(y1,(x-xst)*2+1,1)*bz2)-stamps(2)*ff*
     &                    (vvst_r(y1,(x-xst)*2+1,2)*bbz1
     &                    -vvst_i(y1,(x-xst)*2+1,2)*bbz2)
     &                    )

                     om2i(x,y,z,2)=om2i(x,y,z,2)-fring2(x)*(
     &                    -stamps(1)*ff*
     &                    (vvst_r(y1,(x-xst+1)*2,1)*bz1
     &                    -vvst_i(y1,(x-xst+1)*2,1)*bz2)-stamps(2)*ff*
     &                    (vvst_r(y1,(x-xst+1)*2,2)*bbz1
     &                    -vvst_i(y1,(x-xst+1)*2,2)*bbz2)
     &                    )

                     om2r(x,y,z,3)=om2r(x,y,z,3)-fring1(x)*(
     &                    -stamps(1)*ff*
     &                    (wwst_r(y1,(x-xst)*2+1,1)*bz1
     &                    -wwst_i(y1,(x-xst)*2+1,1)*bz2)-stamps(2)*ff*
     &                    (wwst_r(y1,(x-xst)*2+1,2)*bbz1
     &                    -wwst_i(y1,(x-xst)*2+1,2)*bbz2))

                     om2i(x,y,z,3)=om2i(x,y,z,3)-fring2(x)*(
     &                    -stamps(1)*ff*
     &                    (wwst_r(y1,(x-xst+1)*2,1)*bz1
     &                    -wwst_i(y1,(x-xst+1)*2,1)*bz2)-stamps(2)*ff*
     &                    (wwst_r(y1,(x-xst+1)*2,2)*bbz1
     &                    -wwst_i(y1,(x-xst+1)*2,2)*bbz2))

                  end do
               end do
            end do
         end if
c
c     Fringe forcing for the scalar
c
         if (.not.pert) then
            do ith=1,scalar
               do z=1,nzpc
                  do y=1,min(mby,nyp-yb+1)
                     y1=y+yb-1
                     do x=xst,xen
                        th2r(x,y,z,4+4*(ith-1))=th2r(x,y,z,4+4*(ith-1))-
     &                       fring1(x)*(th2r(x,y,z,1+4*(ith-1))-
     &                       bu1(x,y1,4+(ith-1)))
                        th2i(x,y,z,4+4*(ith-1))=th2i(x,y,z,4+4*(ith-1))-
     &                       fring2(x)*(th2i(x,y,z,1+4*(ith-1))-
     &                       bu2(x,y1,4+(ith-1)))
                     end do
                  end do
               end do
            end do
         else
            do ith=1,scalar
               do z=1,nzpc
                  do y=1,min(mby,nyp-yb+1)
                     y1=y+yb-1
                     do x=xst,xen
                        th2r(x,y,z,4+4*(ith-1))=th2r(x,y,z,4+4*(ith-1))-
     &                       fring1(x)*(th2r(x,y,z,1+4*(ith-1)))
                        th2i(x,y,z,4+4*(ith-1))=th2i(x,y,z,4+4*(ith-1))-
     &                       fring2(x)*(th2i(x,y,z,1+4*(ith-1)))
                     end do
                  end do
               end do
            end do

         end if
      end do

      end subroutine fring
