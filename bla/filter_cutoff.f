c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine filter_cutoff(gu2r,gu2i,ll,npl,nxy,prexn,prezn,wr,wi)
c
c     Spectral cutoff (lowpass) filter at w_c=pi/2
c     filters ll variables of a box in x/z direction
c     input and output is in physical space
c
      implicit none

      include 'par.f'

      integer i,x,y,z,ll,npl,xy,nxy

      real gu2r((nxp/2+1)*mby,nzd,ll)
      real gu2i((nxp/2+1)*mby,nzd,ll)
      real prexn(nx+15),prezn(nz*2+15)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)

      integer xs,xe
      integer zs,ze

      do i=1,ll
c
c     Forward transform
c
         call vrfftf(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &        nx,nzpc*mby,1,nxp/2+1,prexn)
         call vcfftf(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)
c
c     Spectral filter
c
         xs = nx/4+1
         xe = nx/2+1
         do z=1,nz
            do y=1,npl
               do x=xs,xe
                  xy=x+(y-1)*(nxp/2+1)
                  gu2r(xy,z,i) = 0.
                  gu2i(xy,z,i) = 0.
               end do
            end do
         end do
         zs = nz/4+2
         ze = nz-nz/4
         do z=zs,ze
            do y=1,npl
               do x=1,nx/2+1
                  xy=x+(y-1)*(nxp/2+1)
                  gu2r(xy,z,i) = 0.
                  gu2i(xy,z,i) = 0.
               end do
            end do
         end do
c
c     Backtransform and scaling
c
         call vcfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,nz,
     &        nxy,(nxp/2+1)*mby,1,prezn)
         call vrfftb(gu2r(1,1,i),gu2i(1,1,i),wr,wi,
     &        nx,nzpc*mby,1,nxp/2+1,prexn)
         do z=1,nzc
            do y=1,npl
               do x=1,nx/2
                  xy=x+(y-1)*(nxp/2+1)
                  gu2r(xy,z,i) = gu2r(xy,z,i) / real(nx*nz)
                  gu2i(xy,z,i) = gu2i(xy,z,i) / real(nx*nz)
               end do
            end do
         end do
      end do

      end subroutine filter_cutoff
