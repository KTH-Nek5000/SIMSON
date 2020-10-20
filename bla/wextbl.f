c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wextbl(t,vext,cext,re,xl,
     &     xbl,u0low,w0low,longli,dstar,eta,fltype,fend,
     &     my_node,fbla,dybla,nbla)
c
c     Print out the maximum of velocities and vorticities
c     and their coordinates
c
      implicit none

      include 'par.f'
#ifdef MPI
      include 'mpif.h'
#endif

      real t
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real dstar,re,eta(nyp),u0low,w0low,xl
      real xbl
      logical longli
      integer fltype
      real fend
      real ybl,etabl
      real fbla(mbla,7+3*scalar)
      real dybla
      integer nbla

      integer y,i,j
      real vextf(2,7),cextf(2,7,3),vm
c
c     Asymptotic suction boundary layer
c
      real ysuc,usuc,dusuc
      real fext(2,7)
c
c     Functions
c
      real, external :: cubip

c
c     MPI
c
      integer my_node
#ifdef MPI
      integer ip,ierror,nypp,ybp,yb,k
#endif

c
c     Communicate vext and cext
c     Note that the method is very inefficient since every amp(y,i) is
c     sent seperately. Could be improved by defining a strided
c     MPI datatype.
c
#ifdef MPI
      nypp = nyp/nproc+1
      if (my_node.ne.0) then
         do ybp=1,nypp
            yb = (ybp-1)*nproc+my_node+1
            if (yb.le.nyp) then
               do i=1,2
                  do j=1,6
                     call mpi_ssend(vext(yb,i,j),1,mpi_double_precision,
     &                    0,i+10*j+1000*yb,mpi_comm_world,ierror)
                     do k=1,2
                        call mpi_ssend(cext(yb,i,j,k),1,
     &                       mpi_double_precision,
     &                       0,i+10*j+100*k+1000*yb,mpi_comm_world,
     &                       ierror)
                     end do
                  end do
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
                     do i=1,2
                        do j=1,6
                           call mpi_recv(vext(yb,i,j),1,
     &                          mpi_double_precision,ip,i+10*j+1000*yb,
     &                          mpi_comm_world,mpi_status_ignore,ierror)
                           do k=1,2
                              call mpi_recv(cext(yb,i,j,k),1,
     &                             mpi_double_precision,ip,
     &                             i+10*j+100*k+1000*yb,
     &                             mpi_comm_world,mpi_status_ignore,
     &                             ierror)
                           end do
                        end do
                     end do
                  end if
               end do
            end do
         end if
      end if
#endif
c
c     For spatial simulations all maxdata should be refered to the
c     interval fend,xl+fend
c
      if (fltype.ge.4) then
         do i=1,6
            do j=1,2
               do y=1,nyp
                  cext(y,j,i,1)=cext(y,j,i,1)+10.*xl-
     &                 xl*int(10.+(cext(y,j,i,1)-fend*dstar)/xl)
               end do
            end do
         end do
      end if
c
c     This is for Blasius boundary layer
c
      do i=1,7
         vextf(1,i) = 1.E20
         vextf(2,i) =-1.E20
         fext(1,i)  = 1.E20
         fext(2,i)  =-1.E20
      end do
      do y=1,nyp
         if (vextf(1,7).gt.vext(y,1,6)) then
            vextf(1,7)=vext(y,1,6)
            cextf(1,7,1)=cext(y,1,6,1)
            cextf(1,7,2)=eta(y)
            cextf(1,7,3)=cext(y,1,6,2)
         end if
         if (vextf(2,7).lt.vext(y,2,6)) then
            vextf(2,7)=vext(y,2,6)
            cextf(2,7,1)=cext(y,2,6,1)
            cextf(2,7,2)=eta(y)
            cextf(2,7,3)=cext(y,2,6,2)
         end if
c
c     Disturbance amplitude based on the ASBL
c
         ysuc=(eta(y)+1.)/dstar
         usuc=1.-exp(-ysuc)
         dusuc=exp(-ysuc)
         if (fext(1,1).gt.vext(y,1,1)-usuc) then
            fext(1,1)=vext(y,1,1)-usuc
         end if
         if (fext(2,1).lt.vext(y,2,1)-usuc) then
            fext(2,1)=vext(y,2,1)-usuc
         end if
         do i=2,5
            if (fext(1,i).gt.vext(y,1,i)) then
               fext(1,i)=vext(y,1,i)
            end if
            if (fext(2,i).lt.vext(y,2,i)) then
               fext(2,i)=vext(y,2,i)
            end if
         end do
         if (fext(1,6).gt.vext(y,1,6)+dusuc/dstar) then
            fext(1,6)=vext(y,1,6)+dusuc/dstar
         end if
         if (fext(2,6).lt.vext(y,2,6)+dusuc/dstar) then
            fext(2,6)=vext(y,2,6)+dusuc/dstar
         end if

         ybl=1.+eta(y)
         etabl=ybl*sqrt(re/(2.*xbl))

         do i=1,6
c
c     Subtract laminar value
c
            vm=0.
            if (i.eq.1) vm=u0low
            if (i.eq.3) vm=w0low
            if (fltype.eq.3) then
               if (i.eq.1) vm=cubip(etabl,fbla(1,2),dybla,nbla)+u0low
               if (i.eq.6) vm=-cubip(etabl,fbla(1,3),dybla,nbla)*
     &              sqrt(re/(2.*xbl))
            end if
            if (vextf(1,i).gt.vext(y,1,i)-vm) then
               vextf(1,i)=vext(y,1,i)-vm
               cextf(1,i,1)=cext(y,1,i,1)
               cextf(1,i,2)=eta(y)
               cextf(1,i,3)=cext(y,1,i,2)
            end if
            if (vextf(2,i).lt.vext(y,2,i)-vm) then
               vextf(2,i)=vext(y,2,i)-vm
               cextf(2,i,1)=cext(y,2,i,1)
               cextf(2,i,2)=eta(y)
               cextf(2,i,3)=cext(y,2,i,2)
            end if
         end do
      end do
c
c     Write file output "old" output
c
      write(18,*) t/dstar
      write(18,*) fext(1,1), fext(2,1)
      write(18,*) fext(1,2), fext(2,2)
      write(18,*) fext(1,3), fext(2,3)
      write(18,*) fext(1,4)*dstar, fext(2,4)*dstar
      write(18,*) fext(1,5)*dstar, fext(2,5)*dstar
      write(18,*) fext(1,6)*dstar, fext(2,6)*dstar
      do i=1,2
         do j=1,7
            if (j.le.3) then
               write(18,*) vextf(i,j),cextf(i,j,1)/dstar
            else
               write(18,*) vextf(i,j)*dstar,cextf(i,j,1)/dstar
            end if
            write(18,*) (1.+cextf(i,j,2))/dstar,cextf(i,j,3)/dstar
         end do
      end do
      if (longli) then
c
c     y-dependent statistics
c
         do y=1,nyp
            do i=1,2
               do j=1,6
                  if (j.le.3) then
                     write(18,*) vext(y,i,j),cext(y,i,j,1)/dstar
                  else
                     write(18,*) vext(y,i,j)*dstar,cext(y,i,j,1)/dstar
                  end if
                  write(18,*) (1.+eta(y))/dstar,cext(y,i,j,2)/dstar
               end do
            end do
         end do
      end if

#ifdef MPI
      call mpi_barrier(mpi_comm_world,ierror)
#endif

      end subroutine wextbl
