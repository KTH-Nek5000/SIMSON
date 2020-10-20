c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine mhd_poi(ur,ui,zb,zbp,alfa,beta,om3r,om3i,
     &     u3r,u3i,f,fvr,fvi,w3,prey,bis,phr,phi,
     &     q,c,d,e,fmhdr,fmhdi,b0)
c
c     Computes the potential due to the magnetic field imposed
c     by the vector b0.
c     The Lorentz force is added in nonlinbl.
c
c     Note that this subroutine could be positioned in prhs/linearbl to
c     avoid one getxy step.
c
      implicit none
      include 'par.f'
      integer nxz
      parameter (nxz=nx/2*mbz)

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real fmhdr(memnx,memny,memnz,2)
      real fmhdi(memnx,memny,memnz,2)

      real alfa(nxz),beta(nz)
      real om3r(nxz,nyp,3),om3i(nxz,nyp,3)
      real u3r(nxz,nyp,3),u3i(nxz,nyp,3)
      real bbeta(nxz),k2(nxz)
      real fvr(nxz,nyp),fvi(nxz,nyp)
      real f(nxz,nyp,5),bis(nxz,nyp,6)
      real w3(nx/2,mbz,nyp)
      real prey(nyp*2+15)
      real bc(nxz,2,2),fact,bc1
      real phr(nxz,nyp),phi(nxz,nyp)
      real q(nxz,ny+2),c(nxz,ny),d(nxz,ny),e(nxz,ny)
      real temp(nyp),temi(nyp),b0(3)
      integer zb,zbp,x,y,z,xz,i
      integer neumann
      real cc
c
c     zbp is the local index
c     zb is the global index
c
c     streamwise vorticity is in position 4
c     spanwise vorticity is in position 5
c
c     Note: y = -1 is at i = nyp
c           y =  1 is at i = 1

c
c     Boundary condition for potential
c
c     neumann = 0 --> Dirichlet (not used for potential)
c     neumann = 1 --> Neumann (isolated walls)
c     neumann = 2 --> mixed (conducting walls)
c
      neumann = 2
      cc = 0.001
c
c     Get wavenumbers
c
      do z=1,mbz
         do x=1,nx/2
            xz=x+nx/2*(z-1)
            bbeta(xz)=beta(z+zb-1)
            k2(xz)=alfa(xz)*alfa(xz)+bbeta(xz)*bbeta(xz)
         end do
      end do
c
c     Get an xy box of the velocities and vorticities
c
      call getxy(u3r(1,1,1),u3i(1,1,1),zbp,1,ur,ui)
      call getxy(u3r(1,1,2),u3i(1,1,2),zbp,2,ur,ui)
      call getxy(u3r(1,1,3),u3i(1,1,3),zbp,3,ur,ui)

      call getxy(om3r(1,1,1),om3i(1,1,1),zbp,4,ur,ui)
      call getxy(om3r(1,1,3),om3i(1,1,3),zbp,5,ur,ui)
c
c     Construct wall-normal vorticity
c
      do y=1,nyp
         do xz=1,nxz
            om3r(xz,y,2) = -u3i(xz,y,1)*bbeta(xz) +
     &           u3i(xz,y,3)*alfa(xz)
            om3i(xz,y,2) =  u3r(xz,y,1)*bbeta(xz) -
     &           u3r(xz,y,3)*alfa(xz)
         end do
      end do
c
c     Construct right-hand side f = b0_i*omega_i
c     (with normalisation for the Chebyshev transform)
c
      fact = 2./real(nyp-1)
      do y=1,nyp
         do xz=1,nxz
            fvr(xz,y) = fact*(
     &           b0(1)*om3r(xz,y,1) +
     &           b0(2)*om3r(xz,y,2) +
     &           b0(3)*om3r(xz,y,3))
            fvi(xz,y) = fact*(
     &           b0(1)*om3i(xz,y,1) +
     &           b0(2)*om3i(xz,y,2) +
     &           b0(3)*om3i(xz,y,3))
         end do
      end do
c
c     Go into Chebyshev space for the RHS
c     Note that fvr,fvi are mapped onto f(1,1,1) and f(1,1,2)
c
      do i=1,2
         call vchbf(f(1,1,i),w3,nyp,nxz,nxz,1,prey)
      end do

      if (zb.eq.1) then
c
c     For the treatment of the 0/0 mode, store the rhs
c
         do y=1,nyp
            temp(y) = fvr(1,y)
         end do
      end if
c
c     Set the boundary conditions
c
      do i=1,2
         do xz=1,nxz
            bc(xz,1,i) = 0.
            bc(xz,2,i) = 0.
         end do
      end do
c
c     Get the right-handed sides for the Chebyshev tau method
c
      call crhsc(f,bc,nxz,ny,2,nxz,nyp)
c
c     Set up the system and solve for functions (bis contains phr,phi)
c
      call setmac(q,c,d,e,k2,ny,nxz,nxz,neumann,cc)
      call tridc(bis(1,1,1),q,c,d,e,f,nxz,ny,2,nxz,nyp)
c
c     Correction for kx=kz=0 (i.e. Poisson equation)
c
      if (zb.eq.1) then
c
c     Set the boundary condition of the zero/zero mode to zero
c     at the lower wall.
c
c     In principle, the solution from setmac/tridc could be used:
c         phr(1,1) = 0.
c         do y=2,ny-1,2
c            phr(1,1)=phr(1,1)+phr(1,y)-phr(1,y+1)
c         end do
c
c     However, I prefer to actually integrate the rhs directly:
c
c     Set the rows with the boundary conditions to zero (tau method)
c
         temp(ny) = 0.
         temp(ny-1) = 0.
c
c     Integrate twice with zero boundary condition at one boundary.
c     Note that the other Neumann condition is automatically fulfilled
c     in channel flow due to the no-slip conditions.
c
         bc1=0.
         call icheb(temi,temp,bc1,ny,1,1)
         do y=2,ny-1,2
            temi(1)=temi(1)+temi(y)-temi(y+1)
         end do
         call icheb(temp,temi,bc1,ny,1,1)
         do y=2,ny-1,2
            temp(1)=temp(1)+temp(y)-temp(y+1)
         end do
c
c     Overwrite the 0/0 mode with the new values
c
         do y=1,nyp
            phr(1,y) = temp(y)
            phi(1,y) = 0.
         end do

      end if
c
c     Compute wall-normal derivative of phr,phi
c
      call dcheb(bis(1,1,3),phr,nyp,nxz,nxz)
      call dcheb(bis(1,1,4),phi,nyp,nxz,nxz)
c
c     Go back from Chebyshev to physical space
c
      do i=1,4
         call vchbb(bis(1,1,i),w3,nyp,nxz,nxz,1,prey)
      end do
c
c     Now we have:  phr  --> bis(1,1,1)
c                   phi  --> bis(1,1,2)
c                   dphr --> bis(1,1,3)
c                   dphi --> bis(1,1,4)
c
c     Store the potential phi and wall-normal derivative dphi/dy
c
      do i=1,2
         call putxy(bis(1,1,i*2-1),bis(1,1,i*2),zbp,i,fmhdr,fmhdi)
      end do

      if (zb.eq.100000) then
         x = 10
         do y=1,nyp
            write(*,'(4e18.9)')
     &           bis(x,y,1),bis(x,y,2),bis(x,y,3),bis(x,y,4)
         end do
         write(*,*) k2(x)
         write(*,*)
         write(*,*) bis(x,1,3)+cc*k2(x)*bis(x,1,1),
     &        bis(x,1,4)+cc*k2(x)*bis(x,1,2)
         write(*,*) bis(x,ny,3)-cc*k2(x)*bis(x,ny,1),
     &        bis(x,ny,4)-cc*k2(x)*bis(x,ny,2)

         stop
      end if

      end subroutine mhd_poi
