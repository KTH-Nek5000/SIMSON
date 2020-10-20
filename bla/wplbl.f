c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wplbl(npl,nampl,tpl,cpl,
     &     it,b2r,b2i,uxz,uxy,uyz,wr,wi,ur,ui,
     &     mpl,re,xl,zl,t,xs,dstar,fltype,
     &     prexn,prezn,presn,prean,u0low,w0low,alfa,zs,
     &     beta,wplblodd,my_node,realg1,realg2,plxy)

      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
      include 'par.f'
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer z,y,x,yb,i
      character*80 nampl(mpl)
      real re,xl,zl,t,dstar
      integer mpl,tpl(mpl,3),fltype,it,npl
      real cpl(mpl)
      integer ivar,zp,ybp,nypp,myb
      real b2r(nxp/2+1,mby,nzd),b2i(nxp/2+1,mby,nzd)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)
      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real alfa(nx/2*mbz),zs,xs,beta(nz),u0low,w0low
      integer wplblodd(mpl)
      real pi
      parameter (pi = 3.1415926535897932385)
      logical fapp
c
c     Temporary storage arrays (only used on master)
c
      real plxy(nx,nyp/nproc+1)
      real uxy(nx,nyp)
      real uxz(nx,nz)
      real uyz(nyp,nz)
      integer,save :: firsttime=0

c
c     MPI
c
      integer my_node
      integer realg1,realg2
#ifdef MPI
      integer ierror,ip
#endif

      if (nproc.eq.1) then
         nypp = nyp
      else
         nypp = nyp/nproc+1
      end if

      if (my_node.eq.0) then
         write(ios,*) '** write ',npl,' planes (parallel) **',t/dstar
      end if
c
c     Loop over number of planes
c
      do i=1,npl
c
c     On first iteration, compute coordinate and index
c
c         if (it.eq.1) then
         if (firsttime.eq.0) then
            if (tpl(i,1).eq.1) then
c
c     xy plane
c
               cpl(i)   = cpl(i)*dstar
               tpl(i,3) = int(cpl(i)/zl*real(nz)+.5)+nz/2+1
               cpl(i)   = zl/real(nz)*real(tpl(i,3)-nz/2-1)/dstar
               if (tpl(i,3).gt.nz) then
                  write(ioe,*) 'INDEX XY: ',tpl(i,3)
                  call stopnow(545)
               end if
            else if (tpl(i,1).eq.2) then
c
c     xz plane
c
c        The formula for cpl in the Couette and Channel geometries
c       allow for negative y-coordinates


            if (fltype.eq.1.or.fltype.eq.2.or.
     &          fltype.eq.4.or.fltype.eq.5) then

               cpl(i)   = cpl(i)*dstar
               tpl(i,3) = int(acos(cpl(i))/pi*real(nyp-1)+1.5)
               cpl(i)   =(cos(pi*real(tpl(i,3)-1)/real(nyp-1)))/dstar

            else

               cpl(i)   = cpl(i)*dstar-1.
               tpl(i,3) = int(acos(cpl(i))/pi*real(nyp-1)+1.5)
               cpl(i)   =(cos(pi*real(tpl(i,3)-1)/real(nyp-1))+1.)/dstar

            end if


               if (tpl(i,3).gt.nyp) then
                  write(ioe,*) 'INDEX XZ: ',tpl(i,3)
                  call stopnow(545)
               end if


            else if (tpl(i,1).eq.3) then
c
c     yz plane
c
               cpl(i)   = cpl(i)*dstar
               tpl(i,3) = int(cpl(i)/xl*real(nx)+.5)+nx/2+1
               tpl(i,3) = tpl(i,3)-(tpl(i,3)-1)/nx*nx
               wplblodd(i) = mod(tpl(i,3),2)
               tpl(i,3) = (tpl(i,3)+1)/2
               cpl(i)   = int(cpl(i)/xl*real(nx)+.5)*xl/real(nx)
               if (tpl(i,3).gt.nx/2) then
                  write(ioe,*) 'INDEX YZ: ',tpl(i,3)
                  call stopnow(545)
               end if
            else
               call stopnow(454532)
            end if
c
c     Set plane array to zero
c
            uxy = 0.0
c
c     If first iteration, open files and write headers
c     If file exists, do not overwrite it, append the data
c
            if (my_node.eq.0) then
               inquire(file=nampl(i),exist=fapp)
               if(fapp) then
                  write(ios,*) '*** planes appended to file'
                  open(unit=30+i,file=nampl(i),form='unformatted',
     &                 position='append')
               else
                  open(unit=30+i,file=nampl(i),form='unformatted')
                  write(30+i) re,.false.,xl,zl,t,0.
                  write(30+i) nx,nyp,nzc,nfzsym
                  write(30+i) tpl(i,1),tpl(i,2),cpl(i),fltype,dstar
               end if
            end if
         end if

         ivar=tpl(i,2)

         if (tpl(i,1).eq.1) then
c
c     xy plane
c
            zp=tpl(i,3)
c
c     Extract and transform plane i in uxy
c
            do ybp=1,nypp
               yb = (ybp-1)*nproc+my_node+1
               myb= (yb-1)/mby+1
c
c     Get xz plane
c
               if (nproc.eq.1) then
                  call getxz(b2r,b2i,yb,ivar,0,ur,ui)
               else
#ifdef MPI
                  call getpxz(b2r,b2i,yb,ivar,0,ur,ui,
     &                 realg1,realg2,my_node)
#endif
               end if
c
c     Shift and remove moving wall
c
c
c     No shift for Couette flow
c
               if (fltype.ne.2.and.fltype.ne.5) then
                  call xzsh(b2r,b2i,xs,zs,alfa,beta,yb)
                  if (ivar.eq.1) then
                     b2r(1,1,1)=b2r(1,1,1)-u0low
                  end if
                  if (ivar.eq.3) then
                     b2r(1,1,1)=b2r(1,1,1)-w0low
                  end if
               end if
c
c     Transform to physical space
c
               call vcfftb(b2r(1,1,1),b2i(1,1,1),wr,wi,nz,
     &              nx/2,(nxp/2+1)*mby,1,prezn)
               call vrfftb(b2r(1,1,zp),b2i(1,1,zp),wr,wi,
     &              nx,1,1,nxp/2+1,prexn)
c
c     Copy relevant section
c
               do x=1,nx/2
                  plxy(2*x-1,ybp)=b2r(x,1,zp)
                  plxy(2*x,ybp)=b2i(x,1,zp)
               end do
            end do
c
c     Collect the data on processor 0 in uxy
c
            if (my_node.eq.0) then
c
c     Put node 0 into uxy
c
               do y=1,nyp/nproc+1
                  yb=1+(y-1)*nproc
                  if (yb.le.nyp) then
                     do x=1,nx
                        uxy(x,yb)=plxy(x,y)
                     end do
                  end if
               end do
            end if
#ifdef MPI
c
c     Gather information from all processors
c
            if (my_node.gt.0) then
               call mpi_ssend(plxy,nx*(nyp/nproc+1),
     &              mpi_double_precision,
     &              0,my_node+100,mpi_comm_world,ierror)
            else
c
c     Put contribution of other nodes inot uxy
c
               do ip=1,nproc-1
                  call mpi_recv(plxy,nx*(nyp/nproc+1),
     &                 mpi_double_precision,
     &                 ip,ip+100,mpi_comm_world,
     &                 mpi_status_ignore,ierror)
                  do y=1,nyp/nproc+1
                     yb=ip+1+(y-1)*nproc
                     if (yb.le.nyp) then
                        do x=1,nx
                           uxy(x,yb)=plxy(x,y)
                        end do
                     end if
                  end do
               end do
            end if
#endif
            if (my_node.eq.0) then
c
c     Write data
c
               write(30+i) t,0.
               write(30+i) uxy
            end if

         else if (tpl(i,1).eq.2) then
c
c     xz plane
c
            ybp=tpl(i,3)
            yb=(ybp-1)/mby*mby+1+my_node
c
c     get xz plane
c
            if (nproc.eq.1) then
               call getxz(b2r,b2i,yb,ivar,0,ur,ui)
            else
#ifdef MPI
               call getpxz(b2r,b2i,yb,ivar,0,ur,ui,
     &              realg1,realg2,my_node)
#endif
            end if
            if (my_node.eq.0) then
c
c     Shift and remove moving wall
c
               if (fltype.ne.2.and.fltype.ne.5) then
                  call xzsh(b2r,b2i,xs,zs,alfa,beta,yb)
                  if (ivar.eq.1) then
                     b2r(1,1,1)=b2r(1,1,1)-u0low
                  end if
                  if (ivar.eq.3) then
                     b2r(1,1,1)=b2r(1,1,1)-w0low
                  end if
               end if
c
c     Transform to physical space
c
               call vcfftb(b2r(1,1,1),b2i(1,1,1),wr,wi,nz,
     &              nx/2,(nxp/2+1)*mby,1,prezn)
               call vrfftb(b2r(1,1,1),b2i(1,1,1),wr,wi,
     &              nx,nzpc,1,mby*(nxp/2+1),prexn)
c
c     Copy relevant section
c
               do z=1,nz
                  do x=1,nx/2
                     uxz(2*x-1,z)=b2r(x,ybp-yb+1,z)
                     uxz(2*x  ,z)=b2i(x,ybp-yb+1,z)
                  end do
               end do
c
c     Write data
c
               write(30+i) t,0.
               write(30+i) uxz
            end if
         else if (tpl(i,1).eq.3) then
c
c     yz plane
c
            write(ioe,*) 'yz planes not yet implemented with MPI.'
            call stopnow(343454)
         else
            call stopnow(45775)
         end if
      end do

      firsttime = 1


      end subroutine wplbl
