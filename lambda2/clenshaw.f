c ***********************************************************************
C
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine clenshaw(nx,nyp,nypc,nz,ar_in,ai_in,
     &                    ar_out,ai_out,yl,ymin,ymax,dy,fltype,itot)

      implicit none

      integer :: nx,nypc,nz,nyp,itot
      integer :: i,j,jj,k
      real :: aa,bb,yval,yval2,ymax,dy,ymin
      real :: dr(nx/2),di(nx/2),ddr(nx/2),ddi(nx/2),svr(nx/2),svi(nx/2)
      real*8 yl
      real, allocatable :: w3(:,:),prey(:)
      real ar_in(nx/2+1,nyp,nz),ai_in(nx/2+1,nyp,nz)     
      real ar_out(nx/2+1,nypc,nz),ai_out(nx/2+1,nypc,nz)   
      real unif_y(nypc)
      integer :: fltype
c      integer,external :: omp_get_thread_num,omp_get_num_threads

      integer,external :: omp_get_thread_num

      allocate(prey(nyp*2+15))
      call vcosti(nyp,prey,i)
      allocate(w3(nx/2*nz*nyp,itot))
      
!$omp parallel do
      do k=1,nz
          call vchbf(ar_in(1,1,k),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
          call vchbf(ai_in(1,1,k),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
      end do
      deallocate(w3,prey)
c
c     Normalization
c
!$omp parallel do private(i,j,k)
      do k=1,nz
         do j=1,nyp
            do i=1,nx/2+1
               ar_in(i,j,k)=ar_in(i,j,k)*(2./real(nyp-1))
               ai_in(i,j,k)=ai_in(i,j,k)*(2./real(nyp-1))
            end do
         end do
      end do


      if (fltype.eq.1) then
c
c     Channel flow
c
         aa = -1.0
         bb =  1.0
      else
c
c     Boundary layer
c
         aa = 0.0
         bb = yl
         if (ymin < 0.0) then
            unif_y(1) = 0.0
         end if
      end if
      unif_y(1)=ymin
      do j=2,nypc
         unif_y(j)= ymin + real(j-1)/real(nypc-1)*(ymax-ymin)
      end do
      dy = abs(unif_y(2)-unif_y(1))

!$omp parallel do private(i,j,k)
      do k=1,nz
         do j=1,nypc
            do i=1,nx/2+1
               ar_out(i,j,k)=0.
               ai_out(i,j,k)=0.
            end do
         end do
      end do

c
c     Clenshaw recursion
c     


!$omp parallel do private(i,k,jj,yval,yval2,dr,di,ddr,ddi,j,svr,svi)
      do k=1,nz
         do jj=1,nypc
            
            yval = (2.0*unif_y(jj)-aa-bb)/(bb-aa)
            yval2 = 2.0*yval
            
            do i=1,nx/2
               dr(i)  = 0.0
               ddr(i) = 0.0
               di(i)  = 0.0
               ddi(i) = 0.0
            end do
            
            do j=nyp,2,-1
               do i=1,nx/2
                  svr(i) = dr(i)
                  dr(i)  = yval2*dr(i)-ddr(i)+ar_in(i,j,k)
                  ddr(i) = svr(i)
                  
                  svi(i) = di(i)
                  di(i)  = yval2*di(i)-ddi(i)+ai_in(i,j,k)
                  ddi(i) = svi(i)
               end do
            end do

            
            do i=1,nx/2
               ar_out(i,jj,k) = yval*dr(i) - ddr(i) + ar_in(i,1,k)
               ai_out(i,jj,k) = yval*di(i) - ddi(i) + ai_in(i,1,k)
            end do
               
         end do   
      end do
c     
c     The above line is as given by Numerical Recipes
c
c                 ar_out(i,jj,k) = yval*d - dd + 0.5*ar_in(i,1,k)
c
c
c     But the following seems to be the correct formula, 
c     i.e., the first coefficient is not multiplied by half
c
      end subroutine clenshaw
