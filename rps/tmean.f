c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine tmean(umean,tint,nint,tf,tl,itint,pxz,w,uxzs,
     &     npl,sym,nzn)
c
c     Calculate total mean of velocity
c
      implicit none

      include 'par.f'

      integer nzn,npl
      integer ipln,itint(npl),nint
      real tint(npl),tf,tl,umean,sym
      real uxzs(nx,nzn,npl)
      real w(nx,nzn),pxz(nx+2,nz)

      real c,usmean
      integer x,z,it

      umean=0.0
      do ipln=1,nint
         it=itint(ipln)
c
c     Calculate integration weights
c
         call intwgt(c,tint,nint,tf,tl,ipln)
c
c     Add upp statistics
c
         call getxzc(pxz,w,it,sym,npl,uxzs,nzn)
         usmean=0.0
         do z=1,nz
            do x=1,nx
               usmean=usmean+pxz(x,z)
            end do
         end do
         umean=umean+c*usmean/real(nx*nz)
      end do
      return

      end subroutine tmean
