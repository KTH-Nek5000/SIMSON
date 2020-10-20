c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rscale(plane,grid,mxr,mx,my,ivar,ivarth,iscale,
     &     uinf,utau,delta,nux,pdim,pdim2,thdim,thdim2,mxys,mxys2,
     &     mxysth,mxysth2,re,thdelta,tinf,twall)
c
c     Scale variables
c
c     Belongs to  pxyst version 1.0
c     For more info see the pxyst.f file
c
      implicit none

      include 'par.f'

      integer mx,my,ivar,ivarth,iscale,mxr,mxys,mxys2,mxysth,mxysth2
      integer pdim(mxys,2),pdim2(mxys2,2)
      integer thdim(mxysth,3),thdim2(mxysth2,3)
      real re,plane(mx,my),grid(mx,my,2)
      real utau(nx),uinf(nx),delta(nx,3)
      real nux(nx),tinf(nx),twall(nx),thdelta(nx,2)

      integer y,x,pdimv,pdiml,xpr,pdimth
      real c(nx)
      logical plotlog
      character ch

      if (ivar.gt.0) then
         pdimv=pdim(ivar,1)
         pdiml=pdim(ivar,2)
         pdimth=0.
      else
         if (ivar.lt.0) then
            pdimv=pdim2(-ivar,1)
            pdiml=pdim2(-ivar,2)
            pdimth=0.
         else
            if (ivarth.gt.0) then
               pdimv=thdim(ivarth,1)
               pdiml=thdim(ivarth,2)
               pdimth=thdim(ivarth,3)
            else
               if (ivarth.lt.0) then
                  pdimv=thdim2(-ivarth,1)
                  pdiml=thdim2(-ivarth,2)
                  pdimth=thdim2(-ivarth,3)
               end if
            end if
         end if
      end if
      write(*,*) 'Scaling: pdimv,pdiml,pdimth=',pdimv,pdiml,pdimth
c
c     Rescaling dependent variable
c
c     Outer units
c
      if (iscale.eq.1.or.iscale.eq.3.and.ivar.ne.0) then
         do x=1,nx
            c(x)=uinf(x)**pdimv*delta(x,1)**pdiml
         end do
      else
         if (iscale.eq.1.or.iscale.eq.3.and.ivarth.ne.0) then
            do x=1,nx
               c(x)=uinf(x)**pdimv*thdelta(x,1)**pdiml
            end do
         end if
      end if
c
c     Wall units
c
      if (iscale.eq.2.or.iscale.eq.4.or.iscale.eq.5.and.ivar.ne.0) then
         do x=1,nx
            c(x)=utau(x)**pdimv*(1./(re*utau(x)))**pdiml
         end do
      end if
      if (iscale.eq.2.or.iscale.eq.4.or.
     &    iscale.eq.5.and.ivarth.ne.0) then
         do x=1,nx
            c(x)=utau(x)**pdimv*(1./(re*utau(x)))**pdiml
     &           *(nux(x)*abs(twall(x)-tinf(x))/utau(x))**pdimth
         end do
      end if

      do y=1,my
         do x=1,mx
            xpr=mod(10*nx+x+mxr-1,nx)+1
            plane(x,y)=plane(x,y)/c(xpr)
         end do
      end do
c
c     Rescaling y :
c
      if(iscale.eq.3.and.ivar.ne.0) then
         do x=1,nx
            c(x)=delta(x,1)
         end do
      else
         if(iscale.eq.3.and.ivarth.ne.0) then
            do x=1,nx
               c(x)=thdelta(x,1)
            end do
         end if
      end if
      if(iscale.eq.4.or.iscale.eq.5) then
         do x=1,nx
            c(x)=(1./(re*utau(x)))
         end do
      end if
      if(iscale.eq.3.or.iscale.eq.4.or.iscale.eq.5) then
         do y=1,my
            do x=1,mx
               xpr=mod(10*nx+x+mxr-1,nx)+1
               grid(x,y,2)=grid(x,y,2)/c(xpr)
            end do
         end do
      end if
c
c     Plotting logarithmic axis in plus units
c
      plotlog=.false.
      if (iscale.eq.4) then
         write(*,*) 'plot logarithmic in y+ ?'
         read(*,*) ch
         if (ch.eq.'y') then
            plotlog=.true.
         end if
      end if

      if (plotlog) then
         do x=1,mx
            do y=1,my
c               write(*,*) y,grid(x,y,2)
               grid(x,y,2) = log10(grid(x,y,2))
c               write(*,*) y,grid(x,y,2)
            end do
            grid(x,1,2) = grid(x,2,2)
         end do
      end if
c
c     Dividing dependent variable by y+
c
      if (iscale.eq.5) then
         do y=2,my
            do x=1,mx
               plane(x,y)=plane(x,y)/grid(x,y,2)
            end do
         end do
         do x=1,mx
            plane(x,1)=plane(x,2)
         end do
      end if

      end subroutine rscale
