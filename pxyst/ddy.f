c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine ddy(da,a,yl)
c
c     Performs da := ddy (a)
c     the vertical extent of the grid is yl
c     the data is numbered so that y=1 is the lower boundary
c     and y=nyp is the upper (reverse to the normal chebyshev order)
c
      include 'par.f'

      real a(nx,nyp),da(nx,nyp),yl
c
      integer x,y
      real w(nx,nyp)
      real prey(nyp*2+15)
c
c     Initialize transforms and wave numbers
c
      call vcosti(nyp,prey,0)
c
c     Copy the input into the the output array
c     normalize for the transform (2./real(nyp-1))
c     scale for the grid relative to the standard extent  (2./yl)
c     the minus sign accounts for the reverse order
c
      do y=1,nyp
         do x=1,nx
            da(x,y)=a(x,y)*(-2./yl*2./real(nyp-1))
         end do
      end do
c
c     Chebyshev transform
c
      call vchbf(da,w,nyp,nx,nx,1,prey)
c
c     Differentiate in chebyshev space
c
      call rdcheb(da,nyp,nx,nx)
c
c     Return to physical space
c
      call vchbb(da,w,nyp,nx,nx,1,prey)

      return

      end subroutine ddy
