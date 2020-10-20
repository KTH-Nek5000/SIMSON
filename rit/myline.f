c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine myline(graf,gxax,npoint,nplot,sfunc,xmin,xmax,ymin,
     &   ymax,x,z,ix,xp,jvar,sym,my,mys,xl,eta,prex,prez,pxz,w,ur)
c
c     Cuts a line in the y-direction, possibly with subsampling
c
      implicit none

      include 'par.f'

      integer jvar,x,z,my,mys,ix,npoint(20),nplot
      real graf(my,20),gxax(my,20)
      real xp,xmin,xmax,ymin,ymax,sfunc
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15),eta(nyp),xl
      real sym
c
      integer yp,y,x1,i
      real xdmin,c1,c2,fmax
      nplot=(nx-1)/ix+1
c
c     Cut the lines, but wait with scaling and staggering
c
      do yp=1,nyp,mys
         y=(yp-1)/mys+1
c
c     Note how the line is turned upside down to account for reverse numbering
c
         call getxzp(pxz,nyp+1-yp,jvar,ur,sym)
         call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
         call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
         do i=1,nplot
            x1=mod(10*nx+x+ix*(i-1)-1,nx)+1
            graf(y,i)=pxz(x1,z)
            gxax(y,i)=eta(nyp+1-yp)
         end do
      end do
c
c     Decide on scaling x-coordinate=x+c*var
c
      xdmin=0.
      fmax=0.
      do y=1,my
         fmax=max(fmax,abs(graf(y,1)))
         do i=1,nplot-1
            fmax=max(fmax,abs(graf(y,i)))
            xdmin=min(xdmin,graf(y,i+1)-graf(y,i))
         end do
      end do
      if (fmax.eq.0) then
         c1=1.
      else
         if (xdmin.eq.0.) then
            c1=xl/3./fmax
         else
            c1=min(xl/3./fmax,xl/real(nx)*real(ix)/(-1.05*xdmin))
         end if
      end if
c
c     Truncate down to 1,2,4 or 5 * 10**n
c
      c2=10.**int(log10(c1))
      if (log10(c1).lt.0.) c1=c1/10.
      sfunc=5.*c2
      if (c1/c2.lt.5.) sfunc=4.*c2
      if (c1/c2.lt.4.) sfunc=2.*c2
      if (c1/c2.lt.2.) sfunc=c2
c
c     Scale and stagger
c
      xmax=-1.E20
      xmin=1.E20
      do i=1,nplot
         npoint(i)=my
         do y=1,my
            graf(y,i)=xp+real(i-1)*xl/real(nx)*real(ix)+sfunc*
     &           graf(y,i)
            xmax=max(xmax,graf(y,i))
            xmin=min(xmin,graf(y,i))
         end do
      end do
      xmax=xmax+(xmax-xmin)*.05
      xmin=xmin-(xmax-xmin)*.05/1.05
      ymin=eta(nyp)
      ymax=eta(1)

      return

      end subroutine myline
