c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine mypl(ygraph,ygrid,plane,grid,x1,npx,mx,my,mgr)
c
c     Selects data for multiple function of y plots
c
      implicit none

      integer mx,my,mgr,npx,x1(mgr)
      real plane(mx,my),grid(mx,my,2),ygraph(my,mgr),ygrid(my,mgr)
c
      integer i,y
c
      do i=1,npx
         do y=1,my
            ygraph(y,i)=plane(x1(i),y)
            ygrid(y,i)=grid(x1(i),y,2)
         end do
      end do

      return

      end subroutine mypl
