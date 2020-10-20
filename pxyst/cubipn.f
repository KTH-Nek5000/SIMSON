c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cubipn(cip,ipip,xip,nxip,grid,ngrp)
c
c     Calculates coefficients for
c     CUBic InterPolation on Non-equidistant grid
c     This is a slow algorithm, non-vectorized and nxip*ngrp search comparisons
c
c     Input  : xip points to interpolate to
c              nxip number of points in xip
c              grid points to interpolate from
c              ngrp number of grid points
c     Output : cip  cubic interpolation coefficients for each point in xip
c              ipip point first coefficient in cip refers to
c     Usage :  the value of the variable y on grid point xip(k) :
c              is calculated as
c     y(ipip)*cip(k,1)+y(ipip+1)*cip(k,2)+y(ipip+2)*cip(k,3)+y(ipip+3)*cip(k,4)
c
c
      integer nxip,ngrp
      real cip(nxip,4),xip(nxip),grid(ngrp)
      integer ipip(nxip)
c
      integer i
      real dir
      real a(4,4),b(4)
c
c     Find direction of grid
c
      dir=1.
      if (grid(ngrp).lt.grid(1)) dir=-1.
c
c     Find points which surround xip
c
      do i=1,nxip
        if (dir*xip(i).lt.dir*grid(3)) then
          ipip(i)=1
        else
          do j=3,ngrp-2
             if (dir*xip(i).ge.dir*xip(i)) ipip=j-1
          end do
        end if
      end do
c
c     The interpolating polynomial is  :
c     y=b(1)+b(2)*(x-grid(ipip))+b(3)*(x-grid(ipip))**2+b(4)*(x-grid(ipip))**3
c
c     h1=grid(ipip+1)-grid(ipip)
c     h2=grid(ipip+2)-grid(ipip)
c     h3=grid(ipip+3)-grid(ipip)
c     hx=xip-grid(ipip)
c     y1=b(1)
c     y2=b(1)+b(2)*h1+b(3)*h1**2+b(4)*h1**3
c     y3=b(1)+b(2)*h2+b(3)*h2**2+b(4)*h2**3
c     y4=b(1)+b(2)*h3+b(3)*h3**2+b(4)*h3**3
c     yx=b(1)+b(2)*hx+b(3)*hx**2+b(4)*hx**3
c     desired : express  yx in y1 - y4
c
c        (  1  0  0       0      )
c     p= (  1  h1 h1**2   h1**3  )
c        (  1  h2 h2**2   h2**2  )
c        (  1  h3 h3**2   h3**2  )
c
c     y=(y1 y2 y3 y4)^T
c     h=(1 hx hx**2 hx**3)^T
c
c
c     y=pb => b=p^-1 y      b_i=p^-1_ij y_j
c
c     yx = b_i h_i = p^-1_ij y_j b_i = y_j p^-1_ij b_i = y (p^T)^-1 b
c
c     let a=p^T
c
      do i=1,nxip
         do j=1,4
            do k=1,4
               a(j,k)=(grid(ipip(i)+j-1)-grid(ipip(i)))**(k-1)
            end do
         end do
         do k=1,4
            b(j)=(xip(i)-grid(ipip))**(k-1)
         end do
c
c     Solve ax=b
c
      call gauss(a,b)
      end do

      end subroutine cubipn
