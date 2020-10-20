c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine power(n,d,ur,ui,eta,fltype,
     &     prex,prez,pres,prea)
c
c     Adds field n to field 1 with power d
c     If n=1 simply raise to d
c
      implicit none

      include 'par.f'

      integer n,fltype
      real d
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real boxr(nxp/2+1,mby,nzd),boxi(nxp/2+1,mby,nzd)
      real box2r(nxp/2+1,mby,nzd),box2i(nxp/2+1,mby,nzd)
      real eta(nyp)
      real prex(nxp+15),prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)

      integer x,y,z,i,yb,npl
      logical sym
      real wr((nxp/2+1),mby,nzd),wi((nxp/2+1),mby,nzd)
      real um,um1

      if (d.eq.1..and.n.eq.1) return
      do i=1,3
         sym=i.le.2
         do yb=1,nyp,mby
            call getxz(boxr(1,1,1),boxi(1,1,1),yb,i,1,ur,ui)
            npl=min(mby,nyp-yb+1)
            call fft2db(boxr(1,1,1),boxi(1,1,1),sym,npl,
     &           prex,prez,pres,prea,wr,wi)
            if (n.ne.1) then
               call getxz(box2r,box2i,yb,i+3,1,ur,ui)
               call fft2db(box2r(1,1,1),box2i(1,1,1),sym,npl,
     &              prex,prez,pres,prea,wr,wi)
               if (yb.eq.1.and.(fltype.eq.1.or.fltype.eq.4)) then
                  um1=box2r(1,1,1)
               end if
               do y=yb,min(nyp,yb+mby-1)
                  um=eta(y)
                  if (fltype.eq.1.or.fltype.eq.4) um=um1+1.-eta(y)**2
                  if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) um=0.0
                  do z=1,nzpc
                     do x=1,nxp/2
                        boxr(x,y-yb+1,z)=boxr(x,y-yb+1,z)+
     &                       (box2r(x,y-yb+1,z)-um)**d
                        boxi(x,y-yb+1,z)=boxi(x,y-yb+1,z)+
     &                       (box2i(x,y-yb+1,z)-um)**d
                     end do
                  end do
               end do
            else
               do y=yb,min(nyp,yb+mby-1)
                  um=eta(y)
                  if (fltype.eq.1.or.fltype.eq.4) um=um1+1.-eta(y)**2
                  if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) um=0.0
                  do z=1,nzpc
                     do x=1,nxp/2
                        if (boxr(x,y-yb+1,z).lt.0..and.
     &                       boxr(x,y-yb+1,z).gt.-1e-3)
     &                       boxr(x,y-yb+1,z)=0.
                        if (boxi(x,y-yb+1,z).lt.0..and.
     &                       boxi(x,y-yb+1,z).gt.-1e-3)
     &                       boxi(x,y-yb+1,z)=0.
                        boxr(x,y-yb+1,z)=((boxr(x,y-yb+1,z)-um)**d+um)
                        boxi(x,y-yb+1,z)=((boxi(x,y-yb+1,z)-um)**d+um)
                     end do
                  end do
               end do
            end if
            call fft2df(boxr(1,1,1),boxi(1,1,1),sym,npl,
     &           prex,prez,pres,prea,wr,wi)
            do y=yb,min(nyp,yb+mby-1)
               do z=1,nzd
                  do x=1,nxp/2+1
                     boxr(x,y-yb+1,z)=boxr(x,y-yb+1,z)/(nxp*nzp)
                     boxi(x,y-yb+1,z)=boxi(x,y-yb+1,z)/(nxp*nzp)
                  end do
               end do
            end do
            call putxz(boxr,boxi,yb,i,ur,ui)
         end do
      end do

      return

      end subroutine power
