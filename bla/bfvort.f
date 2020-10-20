c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine bfvort(bu1,bu2,bom1,bom2,alfa,prex,prey,prexn,
     &     wbr,wbi,appr,appi,app2r,app2i,w,bu3,bu4,bom3,bom4)
c
c     Computes vorticity associated with base flow
c
      implicit none

      include 'par.f'

      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real bu3(nxp/2+1,nyp,3),bu4(nxp/2+1,nyp,3)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real bom1(nxp/2+1,nyp,3+2*scalar)
      real bom2(nxp/2+1,nyp,3+2*scalar)
      real bom3(nxp/2+1,nyp,3),bom4(nxp/2+1,nyp,3)
      real app1(nyp),app2(nyp)
      real appr(nxp/2+1),appi(nxp/2+1)
      real app2r(nxp/2+1),app2i(nxp/2+1)
      real prey(nyp*2+15),w(nxp/2+1,nyp)
      real alfa(nx/2*mbz),prex(nxp+15),prexn(nx+15)

      integer x,i,y,nn,ith
c
c     y derivatives
c
      nn=0
      do i=1,3,2
         nn=nn+1
         do x=1,nxp/2+1
            do y=1,nyp
               app1(y)=bu1(x,y,i)
               app2(y)=bu2(x,y,i)
            end do
            call vchbf(app1,w,nyp,1,1,1,prey)
            call vchbf(app2,w,nyp,1,1,1,prey)
            call rdcheb(app1,nyp,1,1)
            call rdcheb(app2,nyp,1,1)
            call vchbb(app1,w,nyp,1,1,1,prey)
            call vchbb(app2,w,nyp,1,1,1,prey)
            do y=1,nyp
               bom1(x,y,4-i)=(-1)**nn*app1(y)*(2./real(nyp-1))
               bom2(x,y,4-i)=(-1)**nn*app2(y)*(2./real(nyp-1))
            end do
         end do
      end do
c
c     omy=-dw/dx
c
      do y=1,nyp
         do x=1,nxp/2
            appr(x)=bu1(x,y,3)
            appi(x)=bu2(x,y,3)
         end do
         call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
         do x=1,nxp/2
            app2r(x)=-alfa(x)*appi(x)*(1./real(nxp))
            app2i(x)=alfa(x)*appr(x)*(1./real(nxp))
         end do
         call vrfftb(app2r,app2i,wbr,wbi,nxp,1,1,1,prex)
         do x=1,nxp/2
            bom1(x,y,2)=-app2r(x)
            bom2(x,y,2)=-app2i(x)
         end do
      end do
c
c     omz=omz+dv/dx
c
      do y=1,nyp
         do x=1,nxp/2
            appr(x)=bu1(x,y,2)
            appi(x)=bu2(x,y,2)
         end do
         call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
         do x=1,nxp/2
            app2r(x)=-alfa(x)*appi(x)*(1./real(nxp))
            app2i(x)=alfa(x)*appr(x)*(1./real(nxp))
         end do
         call vrfftb(app2r,app2i,wbr,wbi,nxp,1,1,1,prex)
         do x=1,nxp/2
            bom1(x,y,3)=bom1(x,y,3)+app2r(x)
            bom2(x,y,3)=bom2(x,y,3)+app2i(x)
         end do
      end do
c
c     Do the derivatives of the base flow for the scalars
c
      do ith=1,scalar
c
c     x derivative
c
         do y=1,nyp
            do x=1,nxp/2
               appr(x)=bu1(x,y,3+ith)
               appi(x)=bu2(x,y,3+ith)
            end do
            call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               app2r(x)=-alfa(x)*appi(x)*(1./real(nxp))
               app2i(x)=alfa(x)*appr(x)*(1./real(nxp))
            end do
            call vrfftb(app2r,app2i,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               bom1(x,y,4+2*(ith-1))=app2r(x)
               bom2(x,y,4+2*(ith-1))=app2i(x)
            end do
         end do
c
c     y derivative
c
         do x=1,nxp/2+1
            do y=1,nyp
               app1(y)=bu1(x,y,3+ith)
               app2(y)=bu2(x,y,3+ith)
            end do
            call vchbf(app1,w,nyp,1,1,1,prey)
            call vchbf(app2,w,nyp,1,1,1,prey)
            call rdcheb(app1,nyp,1,1)
            call rdcheb(app2,nyp,1,1)
            call vchbb(app1,w,nyp,1,1,1,prey)
            call vchbb(app2,w,nyp,1,1,1,prey)
            do y=1,nyp
               bom1(x,y,5+2*(ith-1))=app1(y)*(2./real(nyp-1))
               bom2(x,y,5+2*(ith-1))=app2(y)*(2./real(nyp-1))
            end do
         end do


      end do
c
c     Add the periodic point in x
c
      do i=1,3+scalar
         do y=1,nyp
            bom1(nxp/2+1,y,i)=bom1(1,y,i)
            bom2(nxp/2+1,y,i)=bom2(1,y,i)
         end do
      end do

c
c     Find the base flow on the non-dealiasing grid
c     Initialize
c
      do i=1,3
         do y=1,nyp
            do x=1,nxp/2
               bu3(x,y,i)=0.
               bu4(x,y,i)=0.
               bom3(x,y,i)=0.
               bom4(x,y,i)=0.
            end do
         end do
      end do
c
c     Velocities
c
      do i=1,3
         do y=1,nyp
            do x=1,nxp/2
               appr(x)=bu1(x,y,i)
               appi(x)=bu2(x,y,i)
            end do
            call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               appr(x)=appr(x)*(1./real(nxp))
               appi(x)=appi(x)*(1./real(nxp))
            end do
            do x=nx/2+1,nxp/2
               appr(nx/2+1)=0
               appi(nx/2+1)=0
            end do
            call vrfftb(appr,appi,wbr,wbi,nx,1,1,1,prexn)
            do x=1,nx/2
               bu3(x,y,i)=appr(x)
               bu4(x,y,i)=appi(x)
            end do
         end do
      end do
c
c     Vorticities
c
      do i=1,3
         do y=1,nyp
            do x=1,nxp/2
               appr(x)=bom1(x,y,i)
               appi(x)=bom2(x,y,i)
            end do
            call vrfftf(appr,appi,wbr,wbi,nxp,1,1,1,prex)
            do x=1,nxp/2
               appr(x)=appr(x)*(1./real(nxp))
               appi(x)=appi(x)*(1./real(nxp))
            end do
            do x=nx/2+1,nxp/2
               appr(nx/2+1)=0
               appi(nx/2+1)=0
            end do
            call vrfftb(appr,appi,wbr,wbi,nx,1,1,1,prexn)
            do x=1,nx/2
               bom3(x,y,i)=appr(x)
               bom4(x,y,i)=appi(x)
            end do
         end do
      end do

      end subroutine bfvort
