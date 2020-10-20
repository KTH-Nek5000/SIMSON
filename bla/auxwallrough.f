c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine getfxy(boxr,boxi,zb,i,ur,ui)
c
c     To get an xy box from ur and ui.
c
c     size of xy box: (nx/2+1) * mbz * nyp (real and imaginary part)
c
c     The box arrays are on the original (non-dealiased) grid.
c
c     The y-line xi=nx/2+1 is filled with zeros, as done for the z-line
c     xi=nx/2+1 in getxz.
c
      implicit none

      include 'par.f'
c
c     Parameters
c
      integer nxn
      parameter (nxn=nx/2+1)

      integer zb,i
      real boxr(nxn,mbz,nyp),boxi(nxn,mbz,nyp)
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)

      integer x,y,z

      if (mbz.eq.1) then
         do y=1,nyp
            do x=1,nxn-1
               boxr(x,1,y)=ur(x,y,zb,i)
               boxi(x,1,y)=ui(x,y,zb,i)
            end do
            boxr(nxn,1,y)=0.0
            boxi(nxn,1,y)=0.0
         end do
      else
         do z=zb,zb+mbz-1
            do y=1,nyp
               do x=1,nxn-1
                  boxr(x,z-zb+1,y)=ur(x,y,z,i)
                  boxi(x,z-zb+1,y)=ui(x,y,z,i)
               end do
               boxr(nxn,z-zb+1,y)=0.0
               boxi(nxn,z-zb+1,y)=0.0
            end do
         end do
      end if

      end subroutine getfxy
c
c ***********************************************************************
c
      subroutine sfft2df(boxr,boxi,sym,npl,prexn,prezn,presn,prean,
     &     wr,wi)
c
c     To transform npl planes in a box from physical to Fourier space
c     in x & z direction.
c
c     The transforms are performed on the original (non-dealiased) grid
c     ("s" = "small").
c
c     If nfzsym = 1 then select transform symmetry according to flag sym
c
      implicit none

      include 'par.f'
c
c     Further parameters
c
      integer nxn,nzn
      parameter (nxn=nx/2+1,nzn=nz/2+1)

      logical sym
      integer npl
      real boxr(nxn*mby,nz),boxi(nxn*mby,nz)

      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15)
      real prean(nz*3/4+15)
      real wr(nxn*mby,nz),wi(nxn*mby,nz)

      integer nxy
c
c     First real to half complex transform in x-direction
c
      nxy=nxn*npl-1
      call vrfftf(boxr,boxi,wr,wi,nx,nzc*mby,1,nxn,prexn)
c
c     Then complex transform in z-direction
c
      if (nfzsym.eq.0) then
         call vcfftf(boxr,boxi,wr,wi,nz,nxy,nxn*mby,1,prezn)
      else
         if (sym) then
            call vcffts(boxr,boxi,wr,nzn,nxy,nxn*mby,1,presn)
         else
            call vcftaf(boxr(1,2),boxi(1,2),wr,wi,nz/2-1,nxy,
     &           nzn*mby,1,prean)
         end if
      end if

      end subroutine sfft2df
