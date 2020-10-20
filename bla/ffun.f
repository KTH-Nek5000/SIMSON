c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine ffun(fstart,fend,frise,ffall,fmax,
     &     xl,zl,tc,xsc,zsc,
     &     xc1,xc2,fring1,fring2,cphas1,cphas2,sphas1,sphas2,
     &     osnumb,osomega,osalr,osali,osbeta,osamp,osdamp,osmod,my_node)
c
c     Generates the fringe function and computes the phase for
c     Orr-Sommerfeld modes.
c     Computes: xc1, xc2 , fring1, fring2, cphase1, cphase2, sphase1, sphase2
c
      implicit none

      include 'par.f'

      real fstart,fend,frise,ffall,fmax,xl,zl
      real damp1,damp2,phas1,phas2

      real cphas1(osnf,nxp/2,nzd),cphas2(osnf,nxp/2,nzd)
      real sphas1(osnf,nxp/2,nzd),sphas2(osnf,nxp/2,nzd)
      real fring1(nxp/2),fring2(nxp/2)
      real xc1(nxp/2),xc2(nxp/2)

      real fstc,xsc,r,fenc,zc,wt,tc,dnx,zsc,bz
      parameter(r=0.499999999999999999)

      integer xstart,xend,imin,xendc,xst,xen
      integer i,j,x,z

      logical osmod,osdamp
      real osamp
      integer osnumb
      real osbeta(osnf), osomega(osnf)
      real osalr(osnf), osali(osnf)
c
c     Functions
c
      real,external :: step
c
c     MPI
c
      integer my_node

      dnx=xl/real(nxp)

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

         do x=xst,xen
c
c     xc1 is now in channel coordinates
c
            xc1(x)=real(2*x-1-nxp/2-1)*dnx+xsc
            fring1(x)=fmax*(step((xc1(x)-fstc)/frise)-
     &           step((xc1(x)-fenc)/ffall+1.))
            xc2(x)=xc1(x)+dnx
            fring2(x)=fmax*(step((xc2(x)-fstc)/frise)-
     &           step((xc2(x)-fenc)/ffall+1.))
         end do

         if (xstart.gt.xend) then
            do x=xst,xen
               fring1(x)=fring1(x)+fmax
               fring2(x)=fring2(x)+fmax
            end do
         end if

         if (osmod) then

            if (osdamp) then

c -------------------------------------
!$OMP PARALLEL DO PRIVATE(j,wt,bz,x,damp1,damp2,z,zc,phas1,phas2)
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  do x=xst,xen
                     do j=1,osnumb
c
c     If an accurate representation of the modes is necessary,
c     then the full mode with growth/decay should be superimposed:
c
                        damp1 = osamp*exp(-osali(j)*xc1(x))
                        damp2 = osamp*exp(-osali(j)*xc2(x))
c
c     A relative phase between the modes could be added here:
c
c                        if (j.eq.3.or.j.eq.4) then
c                           bz = zc * osbeta(j) + 1.
c                        else
c                           bz = zc * osbeta(j)
c                        end if

                        bz = zc * osbeta(j)
                        wt = tc * osomega(j)
                        phas1 = osalr(j)*xc1(x)+bz-wt
                        phas2 = osalr(j)*xc2(x)+bz-wt
                        cphas1(j,x,z)=cos(phas1)*damp1
                        sphas1(j,x,z)=sin(phas1)*damp1
                        cphas2(j,x,z)=cos(phas2)*damp2
                        sphas2(j,x,z)=sin(phas2)*damp2
                     end do
                  end do
               end do
!$OMP END PARALLEL DO
c -------------------------------------

            else
c
c     However, if the OS/SQ modes contain modes with large decay rates
c     (i.e. large osali) it is inpractical to actually take that
c     large growth in the fringe region (i.e. at positions upstream of
c     the inflow) into account. Then, osali is set to 0:
c
               damp1 = osamp
               damp2 = osamp

c -------------------------------------
!$OMP PARALLEL DO PRIVATE(j,wt,bz,x,z,zc,phas1,phas2)
               do z=1,nzpc
                  zc=zl*real(z-nzp/2-1)/real(nzp)+zsc
                  do x=xst,xen
                     do j=1,osnumb
c
c     A relative phase between the modes could be added here:
c                        bz = zc * osbeta(j) + PHASE
c
                        bz = zc * osbeta(j)
                        wt = tc * osomega(j)
                        phas1 = osalr(j)*xc1(x)+bz-wt
                        phas2 = osalr(j)*xc2(x)+bz-wt
                        cphas1(j,x,z)=cos(phas1)*damp1
                        sphas1(j,x,z)=sin(phas1)*damp1
                        cphas2(j,x,z)=cos(phas2)*damp2
                        sphas2(j,x,z)=sin(phas2)*damp2
                     end do
                  end do
               end do
!$OMP END PARALLEL DO
c -------------------------------------

            end if

         end if
      end do

      end subroutine ffun
