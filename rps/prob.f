c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine prob(plxy,plxz,plyz,upr,y,z,tpl,idxpr,idzpr,
     &     npr,ipln,npl,t,tpr,lsubm,uprm,tint,nint,tf,tl)
c
c     Get probe data from one plane
c
      implicit none

      include 'par.f'

      integer y,z,idxpr,npr,npl,ipln,tpl,idzpr
      real t,upr(npl,100),tpr(npl,100),uprm(100)
      real plxy(nx+1,nyp)
      real plxz(nx+1,nz+1)
      real plyz(nyp,nz+1)
      logical lsubm
      integer nint
      real tf,tl,tint(nint)

      integer i,x,zpr
      real c

      do i=1,npr
         x=mod(1+idxpr*(i-1)-1,nx)+1
         tpr(ipln,i)=t
         if (tpl.eq.1) then
            upr(ipln,i)=plxy(x,nyp+1-y)
         elseif (tpl.eq.2) then
            upr(ipln,i)=plxz(x,z)
         else
            zpr=mod(z+idzpr*(i-1),nz+2)
            upr(ipln,i)=plyz(y,zpr)
         end if
      end do
c
c     Calculate integration weights
c
      if (lsubm) then
         call intwgt(c,tint,nint,tf,tl,ipln)
         do i=1,npr
            if (ipln.eq.1) uprm(i)=0.0
            uprm(i)=uprm(i)+c*upr(ipln,i)
         end do
      end if

      return

      end subroutine prob
