c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine pre_nlin_3d(ur,ui,my_node,prex,pres,prez,prea,
     &     wr,wi,realg1,realg2,u2bf3_r,u2bf3_i,nypp)
c
c     Transform to physical space (vel and vort) and return on the
c     dealiased grid.
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real prex(nxp+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real wr((nxp/2+1)*mby,nzd),wi((nxp/2+1)*mby,nzd)
      real u2bf3_r(nxp/2+1,nzd,nyp/nproc+1,6)
      real u2bf3_i(nxp/2+1,nzd,nyp/nproc+1,6)

      integer ybp,yb,myb,my_node,npl,realg1,realg2,i,nypp
      logical sym

      do ybp=1,nypp
         yb =(ybp-1)*nproc+my_node+1
         myb=(yb-1)/mby+1
         npl=min(mby,nyp-yb+1)

         if (nproc.eq.1) then
            do i=1,6
               call getxz(u2bf3_r(1,1,ybp,i),u2bf3_i(1,1,ybp,i),
     &              yb,i,1,ur,ui)
            end do
         else
#ifdef MPI
            do i=1,6
               call getpxz(u2bf3_r(1,1,ybp,i),u2bf3_i(1,1,ybp,i),yb,
     &              i,1,ur,ui,
     &              realg1,realg2,my_node)
            end do
#endif
         end if

c
c     Backward Fourier transform (Fourier synthesis)
c
         do i=1,3
            sym=i.le.2
            call fft2db(u2bf3_r(1,1,ybp,i),u2bf3_i(1,1,ybp,i),sym,npl,
     &           prex,prez,pres,prea,wr,wi)
            sym=i.eq.3
            call fft2db(u2bf3_r(1,1,ybp,3+i),u2bf3_i(1,1,ybp,3+i),
     &           sym,npl,prex,prez,pres,prea,wr,wi)
         end do

      end do

      end subroutine pre_nlin_3d
