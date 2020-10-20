c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine espec(fesp,esp,upr,npl,npr,tint,nint,tf,tl,
     &     lhanw,nbin,lper)
c
c     Transform probe data to get energy spectra
c
      implicit none

      include 'par.f'

      logical lhanw,lper
      integer nint,npl,npr,nbin
      real upr(npl,100),esp(npl/2,100,2),fesp(npl/2,100)
      real tf,tl,tint(nint)

      integer i,j,k,l
      real domega,c
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     We will do a sloow dft instead of fft since the data
c     is non-equidistant  (if the data is equidistant we could optimize)
c
      domega=2.*pi/(tl-tf)
      do j=1,npr
         do i=1,nint/2
            esp(i,j,1)=0.0
            esp(i,j,2)=0.0
         end do
      end do
c
c     Apply hanging window if desired
c
      if (lhanw) then
         do j=1,nint
            do k=1,npr
               upr(j,k)=upr(j,k)*
     &              .5*(1.+cos(pi*(tint(j)-(tf+tl)*.5)/(tf-tl)*2.))
            end do
         end do
      end if
c
c     First the transform
c
      do j=1,nint
         if (lper) then
            c=1./real(nint)
         else
            call intwgt(c,tint,nint,tf,tl,j)
         end if
         do i=1,nint/2
            do k=1,npr
               esp(i,k,1)=esp(i,k,1)+
     &              c*upr(j,k)*cos(domega*real(i-1)*tint(j))
               esp(i,k,2)=esp(i,k,2)+
     &              c*upr(j,k)*sin(domega*real(i-1)*tint(j))
            end do
         end do
      end do
c
c     Then calculate the energy
c
      do k=1,npr
         esp(1,k,1)=esp(1,k,1)**2
         do i=2,nint/2
            esp(i,k,1)=(esp(i,k,1)**2+esp(i,k,2)**2)*2.
         end do
      end do
c
c     Bin the energy
c
      if (nbin.gt.1) then
         do k=1,npr
            do i=1,nint/2-nbin+1,nbin
               esp(1+(i-1)/nbin,k,1)=esp(i,k,1)
               do l=2,nbin
                  esp(1+(i-1)/nbin,k,1)=esp(1+(i-1)/nbin,k,1)+
     &                                  esp(i+l-1,k,1)
               end do
            end do
         end do
      end if
c
c     Generate frequency axis
c
      do k=1,npr
         do i=1,nint/2
            fesp(i,k)=domega*nbin*real(i-1)
         end do
      end do

      return

      end subroutine espec
