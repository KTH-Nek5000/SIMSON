c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine stagpr(upr,npr,npl,nint,dxpr,xpr,xmoff)
c
c     Stagger the data from multiple probes
c
      implicit none

      real upr(npl,100),xmoff,dxpr,xpr
      integer npr,nint,npl

      integer i,j
      real umin,umax,dxmin
      real urang,c1,c2,stag
      dxmin=0.
      do i=1,npr-1
         do j=1,nint
            dxmin=min(dxmin,upr(j,i+1)-upr(j,i))
         end do
      end do
      urang=0.
      do i=1,npr
         umin=1.e20
         umax=-1.e20
         do j=1,nint
            umin=min(umin,upr(j,i))
            umax=max(umax,upr(j,i))
         end do
         urang=max(urang,umax-umin)
      end do
      stag=max(urang*.5,-dxmin*1.1)
      c1=stag/dxpr
      c2=10.**int(log10(c1))
      if (c1.lt.1.) c2=c2/10.
      xmoff=10.*c2
      if (c1/c2.lt.5.) xmoff=5.*c2
      if (c1/c2.lt.4.) xmoff=4.*c2
      if (c1/c2.lt.2.) xmoff=2.*c2
      do i=1,npr
         do j=1,nint
            upr(j,i)=upr(j,i)+xmoff*(xpr+real(i-1)*dxpr)
         end do
      end do

      return

      end subroutine stagpr
