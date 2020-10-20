c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wampbl(t,it,amp,campw,kx,kz,mwave,nwave,wint,re,
     &     xl,zl,longli,dstar,
     &     rlam,xbl,my_node,eta,fbla,dybla,nbla,spanv,fileurms,fltype,
     &     px,gridy)
c
c     Accumulates and writes the amplitude file
c
c     Function
c
c     amp(y,1)  Streamwise velocity average without (0,0) comp. (squared)
c     amp(y,2)  Normal velocity  average without  (squared)
c     amp(y,3)  Spanwise velocity average without (0,0) (squared)
c     amp(y,4)  Streamwise vorticity average without (0,0) comp. (squared)
c     amp(y,5)  Normal vorticity average (squared)
c     amp(y,6)  Spanwise vorticity average without (0,0) comp.(squared)
c     amp(y,7)  Normal vorticity squared over wavenumber square
c               average, no (0,0)
c     amp(y,8)  Reynolds stress average
c     amp(y,9)  Mean streamwise disturbance velocity
c     amp(y,10) Mean spanwise disturbance velocity
c     amp(y,11) Mean streamwise vorticity component
c     amp(y,12) Mean spanwise vorticity component (to calculate wall shear)
c     amp(y,13-20) empty
c
c     campw(y,i,wave) Complex normal velocity and normal vorticity averages
c                     from selected wavenumbers
c
c     Communicate amp
c     Note that the method is very inefficient since every amp(y,i) is
c     sent seperately. Could be improved by defining a strided
c     MPI datatype.
c     The communication of campw is not yet implemented.
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      integer it,nwave,mwave,kx(nwave),kz(nwave),fltype
      real t,wint(nyp),re,xl,zl,dstar,rlam,xbl,wb
      real amp(nyp,20),spanv
      complex campw(nyp,4,nwave)
      real ybl,etabl
      real eta(nyp),dybla,gridy(nyp)
      real fbla(mbla,7+3*scalar)
      integer nbla
      real px
c
c     Accumulates the amplitudes from the values for each xz-box
c     and writes to logfile and ampfile (unit 15)
c     If longli = .true. then write out y-dependence of statistics
c     first accumulate
c
      logical longli,fileurms
      real sum(100),sumw1,sumw2,sumw3,sumw4,hplus,e0,c,um,omm,dum,wm
      real u0upp,u0low,w0low,tau1,tau2
      integer i,j,y
c
c     Array to store bl thicknesses; delta(1): 99% thickness,
c     delta(2): displacement thickness, delta(3): momentum thickness
c
      real delta(3)
      
c
c     MPI
c
      integer my_node
#ifdef MPI
      integer ip,ierror,nypp,ybp,yb
#endif
c
c     Functions
c
      real,external :: cubip

      real pi
      parameter (pi = 3.1415926535897932385)

#ifdef MPI
      nypp = nyp/nproc+1
      if (my_node.ne.0) then
         do ybp=1,nypp
            yb = (ybp-1)*nproc+my_node+1
            if (yb.le.nyp) then
               do i=1,20
                  call mpi_ssend(amp(yb,i),1,mpi_double_precision,
     &                 0,i+(nproc+30)*yb,mpi_comm_world,ierror)
               end do
            end if
         end do
c
c     Wait here while master is doing the writing
c
         call mpi_barrier(mpi_comm_world,ierror)
         return
      else
         if (nproc.gt.1) then
            do ip=1,nproc-1
               do ybp=1,nypp
                  yb = (ybp-1)*nproc+ip+1
                  if (yb.le.nyp) then
                     do i=1,20
                        call mpi_recv(amp(yb,i),1,mpi_double_precision,
     &                       ip,i+(nproc+30)*yb,mpi_comm_world,
     &                       mpi_status_ignore,ierror)
                     end do
                  end if
               end do
            end do
         end if
      end if
#endif

      u0upp=amp(1,9)
      u0low=amp(nyp,9)
      w0low=amp(nyp,10)

      do j=1,20
         sum(j)=0.
      end do
c
c     This loop makes the amplitude files "compatible" with previous ones
c     by adding in the (0,0) component
c     In addition the turbulence production is calculated
c     this is for Blasius boundary layer
c
      c=-sqrt(re*(rlam+1.)/(2.*xbl))
      do y=1,nyp
         if (
     &        fltype.eq.-2.or.fltype.eq.-1.or.fltype.eq.3.or.
     &        fltype.eq.6.or.fltype.eq.7.or.fltype.eq.8.or.
     &        fltype.eq.9) then
            ybl=1.+eta(y)
            etabl=ybl*sqrt(re/(2.*xbl))
c     When cdev > 0 um has to be recomputed to compare with the (0,0) mode???
            um=cubip(etabl,fbla(1,2),dybla,nbla)+u0low
            wm=spanv*cubip(etabl,fbla(1,5),dybla,nbla)+w0low
            omm=-cubip(etabl,fbla(1,3),dybla,nbla)*sqrt(re/(2.*xbl))
         else if (fltype.eq.1.or.fltype.eq.4) then
            um=1.-eta(y)**2+u0upp
            wm=w0low
            omm=2.*eta(y)
         else if (fltype.eq.2.or.fltype.eq.5) then
            um=(u0upp-u0low)*.5*eta(y)
            wm=w0low
            omm=-(u0upp-u0low)*.5
         else if (fltype.eq.-3) then
c
c     The base flow is not properly defined here!
c
            um = u0low
            wm = w0low
            omm = 0.
         else if (fltype.eq.-20) then
c
c     The base flow is not properly defined here!
c
            um = u0low
            wm = w0low
            omm = 0.
         else
            write(ioe,*) 'Undefined flow type in wampbl subroutine.'
            call stopnow(730113)
         end if

         dum=-omm
c     Comment out amp(:,1) and amp(:,3) to get perturbation velocity only
c         amp(y,1)=amp(y,1)+(amp(y,9)-um)**2
c         amp(y,3)=amp(y,3)+(amp(y,10)-wm)**2
c         amp(y,4)=amp(y,4)+amp(y,11)**2
c         amp(y,6)=amp(y,6)+(amp(y,12)-omm)**2
c         amp(y,7)=amp(y,7)+(amp(y,9)-um)**2+(amp(y,10)-wm)**2
c         amp(y,8)=-dum*amp(y,8)
         amp(y,8) = -amp(y,8)
         amp(y,13) = amp(y,9)
c         amp(y,9)=(amp(y,9)-um)**2
c         amp(y,10)=(amp(y,10)-wm)**2
         amp(y,9 ) = amp(y,9)!**2
         amp(y,10 ) = amp(y,10)!**2
      end do
      do i=1,10
         do y=1,nyp
            sum(i)=sum(i)+amp(y,i)*wint(y)
         end do
      end do
c
c     Compute bulk Reynolds number
c
      do i=13,13
         do y=1,nyp
            sum(i)=sum(i)+amp(y,i)*wint(y)
         end do
      end do
c
c     Compute bl thicknesses
c     delta(1): 99% thickness, delta(2): displacement thickness,
c     delta(3): momentum thickness
c
      do j=1,3
         delta(j)=0.
      end do
      do i=13,13
         do y=1,nyp
            if((delta(1).eq.0.).and.(amp(y,i).le.0.99)) then
               if((amp(y,i).eq.0.99).or.(y.eq.1)) then
                  delta(1)=gridy(y)
               else				
                  delta(1)=gridy(y)+(0.99-amp(y,i))*(gridy(y-1)
     &                 -gridy(y))/(amp(y-1,i)-amp(y,i)) 
               end if
            end if
            delta(2)=delta(2)+(1.-amp(y,i))*wint(y)
            delta(3)=delta(3)+amp(y,i)*(1.-amp(y,i))*wint(y)
         end do
         delta(2)=delta(2)/dstar
         delta(3)=delta(3)/dstar
      end do

c
c     Then scale
c
      do j=1,7
         sum(j)=sqrt(sum(j)*.5)
      end do
      wb=sum(10)*.5
      e0=sqrt((sum(9)+sum(10))*.5)

      tau1 = amp(1,12)
      tau2 = amp(nyp,12)


c
c     For channel and Couette flow
c
      if (fltype.eq.1.or.fltype.eq.2) then
         hplus=sqrt((abs(amp(1,12))+abs(amp(nyp,12)))*re/2.)
      else
c
c     For boundary-layer cases
c
         hplus=sqrt(abs(amp(nyp,12))*re)
      end if
c
c     Write screen output
c
c      write(ios,*) 't',t/dstar,' it ',it
      write(ios,'(A,3(g23.16,tr1))') ' velocity rms            ',
     &                       (sum(j),j=1,3)
      write(ios,'(A,3(g23.16,tr1))') ' vorticity rms           ',
     &                       (sum(j)*dstar,j=4,6)
      write(ios,'(A,4(g23.16,tr1))') ' omy**2/k2, dUuv, e0, h+ ',
     &                       sum(7),sum(8),e0,hplus

c
c     Note that the pressure gradient px might not be
c     correctly computed due to a mismatch of RK substeps etc.
c
      write(ios,'(A,4(g23.16,tr1))') ' px, Reb, tau1/2 ',
     &                       px,sum(13)*re/2.,tau1,tau2
c
c     Write file output
c
c      write(15,*) t/dstar,(sum(j),j=1,3)
c      write(15,*) (sum(j)*dstar,j=4,6),sum(7)
c      write(15,*) sum(8),e0,hplus
      write(15,'(19(e25.16e3,tr1))') t/dstar,(sum(j),j=1,3),
     &     (sum(j)*dstar,j=4,6),sum(7),sum(8),e0,hplus,
     &     px,sum(13)*re/2.,tau1,tau2,(delta(j),j=1,3),wb

      if (longli) then
         open(file='amp.stat-y',position='append',unit=78)
         do i=1,18
            write(78,'(1000e25.16e3)') (amp(y,i),y=1,nyp)
         end do
         close(78)
      end if

      if (fileurms) then
         open(file='energy_urms',position='append',unit=78)
         write(78,'(15e15.7)')t/dstar,(sum(j),j=1,3)
         close(78)
      end if
c
c     Complex amplitudes are saved for selected waves
c
      if (mwave.gt.0) then
#ifdef MPI
c
c     Not implemented for MPI
c
         call stopnow(4343657)
c
c     The following return is to prevent an "internal compiler error"
c     on the XLF compiler
c
         return
#endif
         write(16,*) t/dstar,mwave,nyp,re,2*pi/xl,2*pi/zl,longli
         do i=1,mwave
            c=2.
            if (kx(i).eq.0) c=1.
            sumw1=0.
            sumw2=0.
            sumw3=0.
            sumw4=0.
            do y=1,nyp
               um=0.
               if (kx(i).eq.0.and.kz(i).eq.0) then
                  ybl=1.+eta(y)
                  etabl=ybl*sqrt(re/(2.*xbl))
                  um=cubip(etabl,fbla(1,2),dybla,nbla)+u0low
               end if
               campw(y,1,i)=campw(y,1,i)-um
               sumw1=sumw1+abs(campw(y,1,i))**2*wint(y)*c
               sumw2=sumw2+abs(campw(y,2,i))**2*wint(y)*c
               sumw3=sumw3+abs(campw(y,3,i))**2*wint(y)*c
               sumw4=sumw4+abs(campw(y,4,i))**2*wint(y)*c
            end do
            sumw1=sqrt(sumw1*.5)
            sumw2=sqrt(sumw2*.5)
            sumw3=sqrt(sumw3*.5)
            sumw4=sqrt(sumw4*.5)
            write(16,*) kx(i),kz(i),sumw1,sumw2,sumw3,sumw4
         end do

         if (longli) then
            do i=1,mwave
               if (kx(i).eq.0.and.kz(i).eq.0) then
                  write(16,*) (campw(y,1,i),y=1,nyp)
                  write(16,*) (campw(y,3,i),y=1,nyp)
               else
                  write(16,*) (campw(y,2,i),y=1,nyp)
                  write(16,*) (campw(y,4,i),y=1,nyp)
               end if
            end do
         end if
      end if

#ifdef MPI
      call mpi_barrier(mpi_comm_world,ierror)
#endif

      end subroutine wampbl
