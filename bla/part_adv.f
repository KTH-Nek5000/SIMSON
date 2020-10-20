c************************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine part_adv(eta,an,bn,xl,zl,it,tc,pos,part,pint,pseed,
     &     my_node,fstart,dstar,fltype,taupp,re,bu1,gridy)
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif
      character*80 namn
      integer p,i,ipartstat
      real eta(nyp),an,bn,tc,xl,zl,fstart,dstar,re
      integer it,my_node,fltype
      real part(10,nppart)
      real Pos(3,nppart)
      real uf_1,uf_2,uf
      real vf_1,vf_2,vf
      real wf_1,wf_2,wf
      real y_left,dy
      real appo
      real P_pos(3),P_ux_p,P_uy_p,P_uz_p
      real inv_tau_p,r_part
      integer j_index
      real pint(3,npart,3)
      real taupp(nppart)
      integer, save :: ifileout = 0
      real,external :: pran2
      integer pseed,number_channels
      parameter (number_channels=1)
      real ratio_rho !rho_p/rho_f
      real Re_part
      real bu1(nxp/2+1,nyp,3+scalar)
      real bu1_0(nyp),bv1_0(nyp)
      parameter (ratio_rho=770.)
      real gridy(nyp),grid_Appo(nyp),dy_appo,y_int
      integer x,y
      real u0,v0,cubip
#ifdef MPI
      integer ierror,ip
      integer status1(mpi_status_size)
#endif

c
c     The array Part(npart,10):  3 for velocity, 3 for rhs
c     position, 3 for rhs velocity, 1 tau_p
c
c     Constant
c
      real pi
      parameter (pi = 3.1415926535897932385)




      do y=1,nyp
         bu1_0(y)=bu1(1,y,1)
         bv1_0(y)=bu1(1,y,2)
         grid_appo(y)=acos(gridy(y)*dstar-1.)
      end do

      dy_appo=grid_appo(2)-grid_appo(1)





!$OMP PARALLEL DO PRIVATE(p,appo,j_index,y_left,dy,uf_1,vf_1,wf_1,
!$OMP& uf_2,vf_2,wf_2,inv_tau_p,P_pos,P_ux_p,P_uy_p,P_uz_p,r_part,
!$OMP& uf,vf,wf,Re_part)
      do p=1,nppart


         appo=Pos(2,p)
         j_index = int(acos(appo)/pi*real(nyp-1)+1.)

c            j_index = max(j_index,1)
c            j_index = min(j_index,nyp-1)



         y_left=eta(j_index)
         dy=y_left-eta(j_index+1)
c     Particle radius
c
         r_part=sqrt(4.5*Part(10,p)/(ratio_rho*re))
c
c     Evolve prhs
c
         P_pos(1) = Pos(1,p) + bn * Part(4,p)
         P_pos(2) = Pos(2,p) + bn * Part(5,p)
         P_pos(3) = Pos(3,p) + bn * Part(6,p)
         P_ux_p = Part(1,p) + bn * Part(7,p)
         P_uy_p = Part(2,p) + bn * Part(8,p)
         P_uz_p = Part(3,p) + bn * Part(9,p)

c
c     Interpolation in y
c


         uf_1 = PINT(1,p+my_node*nppart,2-1)
         vf_1 = PINT(2,p+my_node*nppart,2-1)
         wf_1 = PINT(3,p+my_node*nppart,2-1)

         uf_2 = PINT(1,p+my_node*nppart,3-1)
         vf_2 = PINT(2,p+my_node*nppart,3-1)
         wf_2 = PINT(3,p+my_node*nppart,3-1)
c
c     Fluid velocity in particle position
c
         uf = (uf_2-uf_1)*(y_left-Pos(2,p))/dy+uf_1
         vf = (vf_2-vf_1)*(y_left-Pos(2,p))/dy+vf_1
         wf = (wf_2-wf_1)*(y_left-Pos(2,p))/dy+wf_1

c
c     Put in rhs fluid velocity
c


         if (part(10,p).eq.0.) then
c
c     Advect as fluid particle
c
            Part(4,p) = uf
            Part(5,p) = vf
            Part(6,p) = wf
         else
            inv_tau_p=1./Part(10,p)
c
c     Nonlinear Stokes Drag
c
            Re_part=2.*r_part*re*sqrt((uf-Part(1,p))**2+
     &              (vf-Part(2,p))**2+(wf-Part(3,p))**2)


            Part(7,p) = inv_tau_p * (uf - Part(1,p))
     &                  * (1. + 0.15 * Re_part ** 0.687)
            Part(8,p) = inv_tau_p * (vf - Part(2,p))
     &                  * (1. + 0.15 * Re_part ** 0.687)
            Part(9,p) = inv_tau_p * (wf - Part(3,p))
     &                  * (1. + 0.15 * Re_part ** 0.687)


            Part(4,p) = Part(1,p)
            Part(5,p) = Part(2,p)
            Part(6,p) = Part(3,p)
         end if
c---------------------------------------------
c     Save variables for output
c
c     Particle position
c
         pint(1,p,1) = pos(1,p)/dstar
         if (fltype.eq.6) then
            pint(2,p,1) = (pos(2,p)+1)/dstar
         else
            pint(2,p,1) = pos(2,p)/dstar
         end if
         pint(3,p,1) = pos(3,p)/dstar
c
c     Particle velocity
c
         do i=1,3
            pint(i,p,2) = part(i,p)
         end do
c
c     Fluid velocity
c
         pint(1,p,3) = uf
         pint(2,p,3) = vf
         pint(3,p,3) = wf
c----------------------------------------------
c
c     Evolve linear
c


         pos(1,p) = P_pos(1) + an * part(4,p)
         pos(2,p) = P_pos(2) + an * part(5,p)
         pos(3,p) = P_pos(3) + an * part(6,p)

         part(1,p) = P_Ux_p + an * part(7,p)
         part(2,p) = P_Uy_p + an * part(8,p)
         part(3,p) = P_Uz_p + an * part(9,p)
c

c
c     Reposition particles in z (periodic)
c


         if (Pos(3,p).gt.zl/2.)  Pos(3,p)=Pos(3,p)-zl
         if (Pos(3,p).lt.-zl/2.) Pos(3,p)=Pos(3,p)+zl
c
c     Elastic collision at lower wall
c
         if (Pos(2,p).lt.-1.+r_part) then
c
c     Invert y-velocity and correct with particle radius
c
            Part(2,p)=-Part(2,p)
            Pos(2,p) = -2.-Pos(2,p)+2.*r_part


c
c     Deposition
c
c            part(:,p) = 0.
c            pos(2,p) = -1.


         end if
         if (Pos(2,p).gt.1.-r_part) then
            if (fltype.eq.6) then
c
c     Leaving through free stream
c
            Pos(1,p)=0.
            pos(2,p)=((pran2(pseed))*3./2.-1.)
            y_int=acos(pos(2,p))
            u0=cubip(y_int,bu1_0,dy_appo,nyp)
            v0=cubip(y_int,bv1_0,dy_appo,nyp)
            Part(1,p)=u0
            Part(2,p)=v0
            Part(3,p)=0.
            else if (fltype.eq.1.or.fltype.eq.2) then
c
c     Elastic collision at upper wall
c
               Part(2,p)=-Part(2,p)
               Pos(2,p) = 2.-Pos(2,p)-2.*r_part



c
c     Deposition
c
c               part(:,p) = 0.
c               pos(2,p) = 1.

            end if
         end if
c
c       Reposition in x


c
c       For boundary layer particles are setted to x=0
c       with blasius velocity
c
         if((Pos(1,p).gt.xl+fstart).and.fltype.eq.6) then
            Pos(1,p) = 0.
            pos(2,p)=((pran2(pseed))*3./2.-1.)
            y_int=acos(pos(2,p))
            u0=cubip(y_int,bu1_0,dy_appo,nyp)
            v0=cubip(y_int,bv1_0,dy_appo,nyp)
            Part(1,p)=u0
            Part(2,p)=v0
            Part(3,p)=0.
         end if



         if ((fltype.eq.1).or.(fltype.eq.2)) then
            if((number_channels.eq.1).and.
     &        (Pos(1,p).gt.xl)) then
              Pos(1,p) = mod(Pos(1,p),xl)
            end if
            if((number_channels.gt.1).and.
     &        (Pos(1,p).gt.real(number_channels)*xl)) then
              Pos(1,p) = 0.
            end if
         end if
      end do
!$OMP END PARALLEL DO


c
c     Write output for particles
c     Output is written every ipartstat iterations
c
c     We have:
c     pint(:,p,1) particle position
c     pint(:,p,2) particle velocity
c     pint(:,p,3) fluid velocity
c
      ipartstat = 4*50


      if (mod(it-1,ipartstat).eq.0) then
         if (ifileout.eq.0) then
            ! CHECK FOR LOWEST EMPTY FILE if one wants
         end if
         ifileout = ifileout + 1

         if (my_node.eq.0) then
            write(namn,'(a5,i5.5,a5)') 'part-',ifileout,'.stat'
            open(unit=101,file=namn,form='unformatted')
c
c     Write header
c
            write(101) tc/dstar,npart
c
c     Write own data
c
            do p=1,nppart
               write(101)
     &              (pint(i,p,1),i=1,3),(pint(i,p,2),i=1,3),
     &              (pint(i,p,3),i=1,3),part(10,p)
            end do
c
c     Receive from others and write
c
#ifdef MPI
            do ip=1,nproc-1
               call mpi_recv(pint(1,1,1),nppart*3,
     &              mpi_double_precision,
     &              ip,ip+100,mpi_comm_world,status1,ierror)
               call mpi_recv(pint(1,1,2),nppart*3,
     &              mpi_double_precision,
     &              ip,ip+200,mpi_comm_world,status1,ierror)
               call mpi_recv(pint(1,1,3),nppart*3,
     &              mpi_double_precision,
     &              ip,ip+300,mpi_comm_world,status1,ierror)
               call mpi_recv(taupp,nppart,
     &              mpi_double_precision,
     &              ip,ip+400,mpi_comm_world,status1,ierror)

               do p=1,nppart
                  write(101)
     &                 (pint(i,p,1),i=1,3),(pint(i,p,2),i=1,3),
     &                 (pint(i,p,3),i=1,3),taupp(p)
               end do

            end do
#endif
         else
#ifdef MPI
            call mpi_ssend(pint(1,1,1),nppart*3,
     &           mpi_double_precision,0,my_node+100,
     &           mpi_comm_world,ierror)
            call mpi_ssend(pint(1,1,2),nppart*3,
     &           mpi_double_precision,0,my_node+200,
     &           mpi_comm_world,ierror)
            call mpi_ssend(pint(1,1,3),nppart*3,
     &           mpi_double_precision,0,my_node+300,
     &           mpi_comm_world,ierror)
            do p=1,nppart
               taupp(p) = part(10,p)
            end do
            call mpi_ssend(taupp,nppart,
     &           mpi_double_precision,0,my_node+400,
     &           mpi_comm_world,ierror)
#endif
         end if

         if (my_node.eq.0) then
            close(101)
         end if


      end if



      end subroutine part_adv
