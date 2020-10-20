c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine advection(pert,nxy,rot,om2r,om2i,u2r,u2i,
     &     th2r,th2i,lin,bu1,bu2,bom1,bom2,yb,
     &     bf3,u2bf3_r,u2bf3_i,ybp,pr,gr,fltype,dstar)
c
c     Computes the nonlinear terms
c
      implicit none
      include 'par.f'

      integer z,xy,nxy,ith,yb,ybp

      real u2r((nxp/2+1)*mby,nzd,3),u2i((nxp/2+1)*mby,nzd,3)
      real om2r((nxp/2+1)*mby,nzd,3),om2i((nxp/2+1)*mby,nzd,3)
      real th2r((nxp/2+1)*mby,nzd,4*scalar)
      real th2i((nxp/2+1)*mby,nzd,4*scalar)
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bom1(nxp/2+1,nyp,3+2*scalar),bom2(nxp/2+1,nyp,3+2*scalar)
      real u2bf3_r((nxp/2+1)*mby,nzd,nyp/nproc+1,6)
      real u2bf3_i((nxp/2+1)*mby,nzd,nyp/nproc+1,6)
      real h1u,h2u,h1e,h2e,rot
      logical pert,lin,bf3
      real pr(scalar),gr(scalar),fact,chi
      integer fltype
      real pi,dstar
      parameter (pi = 3.1415926535897932385)
c
c     Some parameters for bouyancy boundary layer
c
      chi = 55.*pi/180.
      fact = 1./tan(chi)
c
c     Calculate the full advection term (rotational form)
c     and save it in om2r,om2i
c     in physical space with or without dealiasing
c     Rotational form: H_i = eps_ijk * u_j * omega_k
c     Convective form: N_i = u_j du_i/dx_j
c     Relation:        N_i = .5*d(u_ju_j)/dx_i - H_i
c

      if (.not.pert) then
c
c     Normal nonlinear term for velocities
c
         do z=1,nzpc
            do xy=1,nxy
               h1u = u2r(xy,z,2)*(om2r(xy,z,3)+2.*rot)-
     &              u2r(xy,z,3)*om2r(xy,z,2)
               h2u = u2r(xy,z,3)*om2r(xy,z,1)-
     &              u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)
               om2r(xy,z,3) = u2r(xy,z,1)*om2r(xy,z,2)-
     &              u2r(xy,z,2)*om2r(xy,z,1)
               om2r(xy,z,1) = h1u
               om2r(xy,z,2) = h2u

               h1e = u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)
     &              -u2i(xy,z,3)*om2i(xy,z,2)
               h2e = u2i(xy,z,3)*om2i(xy,z,1)
     &              -u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
               om2i(xy,z,3) = u2i(xy,z,1)*om2i(xy,z,2)
     &              -u2i(xy,z,2)*om2i(xy,z,1)
               om2i(xy,z,1) = h1e
               om2i(xy,z,2) = h2e
            end do
         end do
c
c     Normal nonlinear term for the scalars
c
         do ith=1,scalar
            do z=1,nzpc
               do xy=1,nxy
                  h1u= u2r(xy,z,1) * th2r(xy,z,2+4*(ith-1)) +
     &                 u2r(xy,z,2) * th2r(xy,z,3+4*(ith-1)) +
     &                 u2r(xy,z,3) * th2r(xy,z,4+4*(ith-1))

                  h2u= u2i(xy,z,1) * th2i(xy,z,2+4*(ith-1)) +
     &                 u2i(xy,z,2) * th2i(xy,z,3+4*(ith-1)) +
     &                 u2i(xy,z,3) * th2i(xy,z,4+4*(ith-1))

                  th2r(xy,z,4+4*(ith-1)) = -h1u
                  th2i(xy,z,4+4*(ith-1)) = -h2u
               end do
            end do
         end do
c
c     do forcings related to active scalars
c
         if (abs(fltype).eq.20) then
c
c     Forcing for the buoyancy boundary layer
c
            ith = 1
            do z=1,nzpc
               do xy=1,nxy

                  om2r(xy,z,1) = om2r(xy,z,1) +
     &                 th2r(xy,z,1+4*(ith-1))*gr(ith)/dstar
                  om2i(xy,z,1) = om2i(xy,z,1) +
     &                 th2i(xy,z,1+4*(ith-1))*gr(ith)/dstar

                  om2r(xy,z,2) = om2r(xy,z,2) +
     &                 th2r(xy,z,1+4*(ith-1))*gr(ith)*fact/dstar
                  om2i(xy,z,2) = om2i(xy,z,2) +
     &                 th2i(xy,z,1+4*(ith-1))*gr(ith)*fact/dstar


                  th2r(xy,z,4+4*(ith-1)) = th2r(xy,z,4+4*(ith-1)) -
     &                 4./pr(ith)/gr(ith)*
     &                 (u2r(xy,z,1)+u2r(xy,z,2)*fact)/dstar
                  th2i(xy,z,4+4*(ith-1)) = th2i(xy,z,4+4*(ith-1)) -
     &                 4./pr(ith)/gr(ith)*
     &                 (u2i(xy,z,1)+u2i(xy,z,2)*fact)/dstar

               end do
            end do
         end if
c
c     This is for an active scalar
c
         do ith=1,scalar
            if (gr(ith).ne.0.) then
               do z=1,nzpc
                  do xy=1,nxy
                     om2r(xy,z,2) = om2r(xy,z,2) +
     &                    th2r(xy,z,1+4*(ith-1))*gr(ith)
                     om2i(xy,z,2) = om2i(xy,z,2) +
     &                    th2i(xy,z,1+4*(ith-1))*gr(ith)
                  end do
               end do
            end if
         end do

      else

         if(.not.bf3) then
            if (lin) then
c
c     Perturbation, linear, 2D base flow
c
               do z=1,nzpc
                  do xy=1,nxy
                     h1u= u2r(xy,z,2)*(bom1(xy,yb,3)+2.*rot)-
     &                    u2r(xy,z,3)*bom1(xy,yb,2)+
     &                    bu1(xy,yb,2)*(om2r(xy,z,3)+2.*rot)-
     &                    bu1(xy,yb,3)*om2r(xy,z,2)

                     h2u= u2r(xy,z,3)*bom1(xy,yb,1)-
     &                    u2r(xy,z,1)*(bom1(xy,yb,3)+2.*rot)+
     &                    bu1(xy,yb,3)*om2r(xy,z,1)-
     &                    bu1(xy,yb,1)*(om2r(xy,z,3)+2.*rot)

                     om2r(xy,z,3)=u2r(xy,z,1)*bom1(xy,yb,2)-
     &                    u2r(xy,z,2)*bom1(xy,yb,1)+
     &                    bu1(xy,yb,1)*om2r(xy,z,2)-
     &                    bu1(xy,yb,2)*om2r(xy,z,1)

                     om2r(xy,z,1)=h1u
                     om2r(xy,z,2)=h2u

                     h1e= u2i(xy,z,2)*(Bom2(xy,yb,3)+2.*rot)-
     &                    u2i(xy,z,3)*Bom2(xy,yb,2)+
     &                    bu2(xy,yb,2)*(om2i(xy,z,3)+2.*rot)-
     &                    bu2(xy,yb,3)*om2i(xy,z,2)

                     h2e= u2i(xy,z,3)*Bom2(xy,yb,1)-
     &                    u2i(xy,z,1)*(Bom2(xy,yb,3)+2.*rot)+
     &                    bu2(xy,yb,3)*om2i(xy,z,1)-
     &                    bu2(xy,yb,1)*(om2i(xy,z,3)+2.*rot)

                     om2i(xy,z,3)=u2i(xy,z,1)*Bom2(xy,yb,2)-
     &                    u2i(xy,z,2)*Bom2(xy,yb,1)+
     &                    bu2(xy,yb,1)*om2i(xy,z,2)-
     &                    bu2(xy,yb,2)*om2i(xy,z,1)

                     om2i(xy,z,1) = h1e
                     om2i(xy,z,2) = h2e
                  end do
               end do
c
c     Perturbation, linear, 2D base flow, for the scalar
c
               do ith=1,scalar
                  do z=1,nzpc
                     do xy=1,nxy

                        h1u= bu1(xy,yb,1) * th2r(xy,z,2+4*(ith-1)) +
     &                       bu1(xy,yb,2) * th2r(xy,z,3+4*(ith-1)) +
     &                       bu1(xy,yb,3) * th2r(xy,z,4+4*(ith-1)) +
     &                       u2r(xy,z,1) * bom1(xy,yb,4+2*(ith-1)) +
     &                       u2r(xy,z,2) * bom1(xy,yb,5+2*(ith-1))

                        h2u= bu2(xy,yb,1) * th2i(xy,z,2+4*(ith-1)) +
     &                       bu2(xy,yb,2) * th2i(xy,z,3+4*(ith-1)) +
     &                       bu2(xy,yb,3) * th2i(xy,z,4+4*(ith-1)) +
     &                       u2i(xy,z,1) * bom2(xy,yb,4+2*(ith-1)) +
     &                       u2i(xy,z,2) * bom2(xy,yb,5+2*(ith-1))

                        th2r(xy,z,4+4*(ith-1)) = -h1u
                        th2i(xy,z,4+4*(ith-1)) = -h2u

                     end do
                  end do
               end do


            else
c
c     Perturbation, nonlinear, 2D base flow
c
               do z=1,nzpc
                  do xy=1,nxy
                     h1u= u2r(xy,z,2)*(bom1(xy,yb,3)+2.*rot)-
     &                    u2r(xy,z,3)* bom1(xy,yb,2)+
     &                    bu1(xy,yb,2)*(om2r(xy,z,3)+2.*rot)-
     &                    bu1(xy,yb,3)*om2r(xy,z,2)+
     &                    u2r(xy,z,2)*(om2r(xy,z,3)+2.*rot)-
     &                    u2r(xy,z,3)*om2r(xy,z,2)

                     h2u= u2r(xy,z,3)*bom1(xy,yb,1)-
     &                    u2r(xy,z,1)*(bom1(xy,yb,3)+2.*rot)+
     &                    bu1(xy,yb,3)*om2r(xy,z,1)-
     &                    bu1(xy,yb,1)*(om2r(xy,z,3)+2.*rot)+
     &                    u2r(xy,z,3)*om2r(xy,z,1)-
     &                    u2r(xy,z,1)*(om2r(xy,z,3)+2.*rot)

                     om2r(xy,z,3)=u2r(xy,z,1)*bom1(xy,yb,2)-
     &                    u2r(xy,z,2)*bom1(xy,yb,1)+
     &                    bu1(xy,yb,1)*om2r(xy,z,2)-
     &                    bu1(xy,yb,2)*om2r(xy,z,1)+
     &                    u2r(xy,z,1)*om2r(xy,z,2)-
     &                    u2r(xy,z,2)*om2r(xy,z,1)

                     om2r(xy,z,1)=h1u
                     om2r(xy,z,2)=h2u

                     h1e= u2i(xy,z,2)*(bom2(xy,yb,3)+2.*rot)-
     &                    u2i(xy,z,3)*bom2(xy,yb,2)+
     &                    bu2(xy,yb,2)*(om2i(xy,z,3)+2.*rot)-
     &                    bu2(xy,yb,3)*om2i(xy,z,2)+
     &                    u2i(xy,z,2)*(om2i(xy,z,3)+2.*rot)-
     &                    u2i(xy,z,3)*om2i(xy,z,2)

                     h2e= u2i(xy,z,3)*bom2(xy,yb,1)-
     &                    u2i(xy,z,1)*(bom2(xy,yb,3)+2.*rot)+
     &                    bu2(xy,yb,3)*om2i(xy,z,1)-
     &                    bu2(xy,yb,1)*(om2i(xy,z,3)+2.*rot)+
     &                    u2i(xy,z,3)*om2i(xy,z,1)-
     &                    u2i(xy,z,1)*(om2i(xy,z,3)+2.*rot)
                     om2i(xy,z,3)=u2i(xy,z,1)*bom2(xy,yb,2)-
     &                    u2i(xy,z,2)*bom2(xy,yb,1)+
     &                    bu2(xy,yb,1)*om2i(xy,z,2)-
     &                    bu2(xy,yb,2)*om2i(xy,z,1)+
     &                    u2i(xy,z,1)*om2i(xy,z,2)-
     &                    u2i(xy,z,2)*om2i(xy,z,1)

                     om2i(xy,z,1) = h1e
                     om2i(xy,z,2) = h2e
                  end do
               end do
c
c     Perturbation, nonlinear, 2D base flow, for the scalar
c
               do ith=1,scalar
                  do z=1,nzpc
                     do xy=1,nxy

                        h1u= u2r(xy,z,1) * th2r(xy,z,2+4*(ith-1)) +
     &                       u2r(xy,z,2) * th2r(xy,z,3+4*(ith-1)) +
     &                       u2r(xy,z,3) * th2r(xy,z,4+4*(ith-1)) +
     &                       bu1(xy,yb,1) * th2r(xy,z,2+4*(ith-1)) +
     &                       bu1(xy,yb,2) * th2r(xy,z,3+4*(ith-1)) +
     &                       bu1(xy,yb,3) * th2r(xy,z,4+4*(ith-1)) +
     &                       u2r(xy,z,1) * bom1(xy,yb,4+2*(ith-1)) +
     &                       u2r(xy,z,2) * bom1(xy,yb,5+2*(ith-1))

                        h2u= u2i(xy,z,1) * th2i(xy,z,2+4*(ith-1)) +
     &                       u2i(xy,z,2) * th2i(xy,z,3+4*(ith-1)) +
     &                       u2i(xy,z,3) * th2i(xy,z,4+4*(ith-1)) +
     &                       bu2(xy,yb,1) * th2i(xy,z,2+4*(ith-1)) +
     &                       bu2(xy,yb,2) * th2i(xy,z,3+4*(ith-1)) +
     &                       bu2(xy,yb,3) * th2i(xy,z,4+4*(ith-1)) +
     &                       u2i(xy,z,1) * bom2(xy,yb,4+2*(ith-1)) +
     &                       u2i(xy,z,2) * bom2(xy,yb,5+2*(ith-1))

                        th2r(xy,z,4+4*(ith-1)) = -h1u
                        th2i(xy,z,4+4*(ith-1)) = -h2u

                     end do
                  end do
               end do

            end if



c
c     do forcings related to active scalars, perturbation form
c
            if (abs(fltype).eq.20) then
c
c     Forcing for the buoyancy boundary layer
c
               ith = 1
               do z=1,nzpc
                  do xy=1,nxy

                     om2r(xy,z,1) = om2r(xy,z,1) +
     &                    th2r(xy,z,1+4*(ith-1))*gr(ith)/dstar
                     om2i(xy,z,1) = om2i(xy,z,1) +
     &                    th2i(xy,z,1+4*(ith-1))*gr(ith)/dstar

                     om2r(xy,z,2) = om2r(xy,z,2) +
     &                    th2r(xy,z,1+4*(ith-1))*gr(ith)*fact/dstar
                     om2i(xy,z,2) = om2i(xy,z,2) +
     &                    th2i(xy,z,1+4*(ith-1))*gr(ith)*fact/dstar

                     th2r(xy,z,4+4*(ith-1)) = th2r(xy,z,4+4*(ith-1))-
     &                    4./pr(ith)/gr(ith)*
     &                    (u2r(xy,z,1)+u2r(xy,z,2)*fact)/dstar
                     th2i(xy,z,4+4*(ith-1)) = th2i(xy,z,4+4*(ith-1))-
     &                    4./pr(ith)/gr(ith)*
     &                    (u2i(xy,z,1)+u2i(xy,z,2)*fact)/dstar

                  end do
               end do
            else if (scalar.gt.0.and.gr(1).ne.0) then
               call stopnow(999)
c
c     This is for an active scalar
c
c     do ith=1,scalar
c     do z=1,nzpc
c     do xy=1,nxy
c     if (gr(ith).ne.0.) then
c     om2r(xy,z,2) = om2r(xy,z,2) +
c     &                          th2r(xy,z,1+4*(ith-1))*gr(ith)
c     om2i(xy,z,2) = om2i(xy,z,2) +
c     &                          th2i(xy,z,1+4*(ith-1))*gr(ith)
c     end if
c     end do
c     end do
c     end do
            end if

         else

            if (scalar.gt.0) then
               write(*,*) 'scalar and 3D base flow not implemented'
               call stopnow(435452)
            end if


            if (lin) then
c
c     Perturbation, linearized, 3D base flow
c
               do z=1,nzpc
                  do xy=1,nxy
                     h1u=     u2r(xy,z,2)*(u2bf3_r(xy,z,ybp,6)
     &                    +2.*rot)-
     &                    u2r(xy,z,3)*u2bf3_r(xy,z,ybp,5)+
     &                    u2bf3_r(xy,z,ybp,2)*(om2r(xy,z,3))-
     &                    u2bf3_r(xy,z,ybp,3)*om2r(xy,z,2)

                     h2u=     u2r(xy,z,3)*u2bf3_r(xy,z,ybp,4)-
     &                    u2r(xy,z,1)*(u2bf3_r(xy,z,ybp,6)+2.*rot)+
     &                    u2bf3_r(xy,z,ybp,3)*om2r(xy,z,1)-
     &                    u2bf3_r(xy,z,ybp,1)*(om2r(xy,z,3))

                     om2r(xy,z,3)=u2r(xy,z,1)*u2bf3_r(xy,z,ybp,5)-
     &                    u2r(xy,z,2)*u2bf3_r(xy,z,ybp,4)+
     &                    u2bf3_r(xy,z,ybp,1)*om2r(xy,z,2)-
     &                    u2bf3_r(xy,z,ybp,2)*om2r(xy,z,1)


                     om2r(xy,z,1)=h1u
                     om2r(xy,z,2)=h2u

                     h1e=     u2i(xy,z,2)*(u2bf3_i(xy,z,ybp,6)
     &                    +2.*rot)-
     &                    u2i(xy,z,3)*u2bf3_i(xy,z,ybp,5)+
     &                    u2bf3_i(xy,z,ybp,2)*(om2i(xy,z,3))-
     &                    u2bf3_i(xy,z,ybp,3)*om2i(xy,z,2)

                     h2e=     u2i(xy,z,3)*u2bf3_i(xy,z,ybp,4)-
     &                    u2i(xy,z,1)*(u2bf3_i(xy,z,ybp,6)+2.*rot)+
     &                    u2bf3_i(xy,z,ybp,3)*om2i(xy,z,1)-
     &                    u2bf3_i(xy,z,ybp,1)*(om2i(xy,z,3))

                     om2i(xy,z,3)=u2i(xy,z,1)*u2bf3_i(xy,z,ybp,5)-
     &                    u2i(xy,z,2)*u2bf3_i(xy,z,ybp,4)+
     &                    u2bf3_i(xy,z,ybp,1)*om2i(xy,z,2)-
     &                    u2bf3_i(xy,z,ybp,2)*om2i(xy,z,1)

                     om2i(xy,z,1) = h1e
                     om2i(xy,z,2) = h2e
                  end do
               end do
            else
c
c     Perturbation, nonlinear, 3D base flow
c
               do z=1,nzpc
                  do xy=1,nxy
                     h1u=     u2r(xy,z,2)*(u2bf3_r(xy,z,ybp,6)
     &                    +2.*rot)-
     &                    u2r(xy,z,3)*u2bf3_r(xy,z,ybp,5)+
     &                    u2bf3_r(xy,z,ybp,2)*(om2r(xy,z,3))-
     &                    u2bf3_r(xy,z,ybp,3)*om2r(xy,z,2)+
     &                    u2r(xy,z,2)*(om2r(xy,z,3))-
     &                    u2r(xy,z,3)*om2r(xy,z,2)

                     h2u=     u2r(xy,z,3)*u2bf3_r(xy,z,ybp,4)-
     &                    u2r(xy,z,1)*(u2bf3_r(xy,z,ybp,6)+2.*rot)+
     &                    u2bf3_r(xy,z,ybp,3)*om2r(xy,z,1)-
     &                    u2bf3_r(xy,z,ybp,1)*(om2r(xy,z,3))+
     &                    u2r(xy,z,3)*om2r(xy,z,1)-
     &                    u2r(xy,z,1)*(om2r(xy,z,3))

                     om2r(xy,z,3)=u2r(xy,z,1)*u2bf3_r(xy,z,ybp,5)-
     &                    u2r(xy,z,2)*u2bf3_r(xy,z,ybp,4)+
     &                    u2bf3_r(xy,z,ybp,1)*om2r(xy,z,2)-
     &                    u2bf3_r(xy,z,ybp,2)*om2r(xy,z,1)+
     &                    u2r(xy,z,1)*om2r(xy,z,2)-
     &                    u2r(xy,z,2)*om2r(xy,z,1)

                     om2r(xy,z,1)=h1u
                     om2r(xy,z,2)=h2u

                     h1e=     u2i(xy,z,2)*(u2bf3_i(xy,z,ybp,6)
     &                    +2.*rot)-
     &                    u2i(xy,z,3)*u2bf3_i(xy,z,ybp,5)+
     &                    u2bf3_i(xy,z,ybp,2)*(om2i(xy,z,3))-
     &                    u2bf3_i(xy,z,ybp,3)*om2i(xy,z,2)+
     &                    u2i(xy,z,2)*(om2i(xy,z,3))-
     &                    u2i(xy,z,3)*om2i(xy,z,2)

                     h2e=     u2i(xy,z,3)*u2bf3_i(xy,z,ybp,4)-
     &                    u2i(xy,z,1)*(u2bf3_i(xy,z,ybp,6)+2.*rot)+
     &                    u2bf3_i(xy,z,ybp,3)*om2i(xy,z,1)-
     &                    u2bf3_i(xy,z,ybp,1)*(om2i(xy,z,3))+
     &                    u2i(xy,z,3)*om2i(xy,z,1)-
     &                    u2i(xy,z,1)*(om2i(xy,z,3))
                     om2i(xy,z,3)=u2i(xy,z,1)*u2bf3_i(xy,z,ybp,5)-
     &                    u2i(xy,z,2)*u2bf3_i(xy,z,ybp,4)+
     &                    u2bf3_i(xy,z,ybp,1)*om2i(xy,z,2)-
     &                    u2bf3_i(xy,z,ybp,2)*om2i(xy,z,1)+
     &                    u2i(xy,z,1)*om2i(xy,z,2)-
     &                    u2i(xy,z,2)*om2i(xy,z,1)

                     om2i(xy,z,1) = h1e
                     om2i(xy,z,2) = h2e
                  end do
               end do
            end if
         end if
      end if

      end subroutine advection
