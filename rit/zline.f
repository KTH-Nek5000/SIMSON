c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine zline(graf,gxax,imode,iopt,xmin,xmax,ymin,ymax,x,y,
     &     jvar,sym,mz,mzs,zl,prex,prez,pxz,w,ur,inunit,utunit)
c
      implicit none

      include 'par.f'

      integer jvar,x,z,y,mz,mzs,iopt,inunit,utunit,imode
      integer xpr
      real graf(mz),gxax(mz)
      real xmin,xmax,ymin,ymax
      real pxz(nx+2,nz),w(nx+2,nz)
      complex ur(memnx,memny,memnz,4)
      real prex(nx+15),prez(nz*2+15)
      real sym,zl
      character*1 ans
c
c     Cut the line
c
      call getxzp(pxz,y,jvar,ur,sym)
      call vcfftb(pxz,pxz(2,1),w,w(2,1),nz,nx/2,nx+2,2,prez)
      call vrfftb(pxz,pxz(2,1),w,w(2,1),nx,nz,2,nx+2,prex)
      xpr=mod(10*nx+x-1,nx)+1
      do z=1,nz,mzs
         graf(z)=pxz(xpr,z)
         gxax(z)=zl/real(nz)*real(z-nz/2-1)
      end do
c
c     Extend periodically in the z-direction
c
      graf(nz/mzs+1)=graf(1)
      gxax(nz/mzs+1)=zl/real(nz)*real(nz+mzs-nz/2-1)
      ymin=0.
      ymax=0.
      do z=1,nz/mzs+1
         ymin=min(ymin,graf(z))
         ymax=max(ymax,graf(z))
      end do
      ymin=ymin*1.2
      ymax=ymax*1.2
      xmin=gxax(1)
      xmax=gxax(nz/mzs+1)

      write(*,*) 'do you want to get file'
      read(inunit,9000) ans
      if (ans.eq.'Y'.or.ans.eq.'y') then
         write(55,*)'x=',x,'y=',y
         do z=1,nz,mzs
            write(55,*) gxax(z),graf(z)
         end do
      end if
      if (iopt.eq.1) then
         write(*,*) 'Do you want to change ymin and ymax (y/n)'
         read(inunit,9000) ans
         if (imode.ge.2) write(utunit,9000) ans
 9000    format(a1)
         if (ans.eq.'Y'.or.ans.eq.'y') then
            if (imode.eq.3) then
               write(*,*)' Give ymin and ymax '
            else
               write(*,*)' Give ymin and ymax ( now ',ymin,ymax,')'
            end if
            read(inunit,*) ymin,ymax
            if (imode.ge.2) write(utunit,*) ymin,ymax
         end if
      end if

      return

      end
