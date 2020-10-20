c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine foustr(xxc,maxmod,xr,xl,heada,headb,re,mbox,ihard)
c
c     Plots the curve of the maximum over y
c     of the choosen variable and growth rate
c
c
c
      implicit none

      include 'par.f'

      character*1 ans
      integer mbox,ihard
      character*80 heada,headb,headc
      real xr,xl,re
      real xbr(nx),ybr(nx),y2(nx+nyp)
      real maxmod(0:nx,3),xxc(0:nx,20)
      real ux(nx),umax(nx),uxx(nx)
      integer i,j,x,nbr,npoint(1)
      integer deg
      real ampx(4*nx),coeff(10),xx(4*nx,1),yy(4*nx,1)
      real start,end,ys,ye,amp(4*nx),damp(4*nx)


 1    do x=1,nx
         ux(x)=float(x-nx/2-1)/float(nx)*xl+xr
         uxx(x)=ux(x)
c         write(51,*)ux(x),xxc(x-1,1)
c         ux(x)=1.720787657571*sqrt(re*(ux(x)+re/1.720787657571**2))
         umax(x)=maxmod(x-1,1)
c         write(51,*)ux(x),xxc(x-1,1)
      end do

      write(*,9200) uxx(1),uxx(nx)
 9200 format('give start and end for x (',f7.2,'-',f7.2,')')
      read(*,*) start,end
c     start=1.720787657571*sqrt(re*(start+re/1.720787657571**2))
c     end=1.720787657571*sqrt(re*(end+re/1.720787657571**2))
      if (start.lt.ux(1).or.start.gt.ux(nx)) start=ux(1)
      if (end.gt.ux(nx).or.end.lt.ux(1)) end=ux(nx)
      ys=1
      ye=1
      mbox=1

      write(*,9300) nx,4*nx
 9300 format('give number of re1 points (0 for default ',i5,
     &     ' max ',i5,')')
      read(*,*) npoint(1)
      if (npoint(1).eq.0) npoint(1)=nx
      do i=1,npoint(1)
         ampx(i)=start+i*(end-start)/npoint(1)
      end do

      write(*,*) ' degree of polynomial'
      read(*,*) deg

      if (deg.gt.0) then
         j=1
         do i=1,nx
            if (ux(i).ge.start.and.ux(i).le.end) then
               ux(j)=ux(i)
               umax(j)=umax(i)
               j=j+1
            end if
         end do
        call curvefit(ux,umax,nx,j-1,deg,ampx,amp,damp,npoint(1),coeff)
        write(headc,4000)deg
 4000   format('Degree of polynomial =',i4)
      else
         headc='Spline'
         call spline(ux,umax,y2,nx,xbr,ybr,nbr)
         call splint(ux,umax,y2,nx,ampx,amp,damp,npoint(1))
      end if

      if (amp(npoint(1)).lt.0.) then
         write(*,*) '    '
         write(*,*) '**  Try with lower degre of polynomial  **'
         write(*,*) '    '
         write(*,*) '    '
         write(*,*) '    '
         return
      end if
      do i=1,npoint(1)
         xx(i,1)=0
         yy(i,1)=0
      end do

      do i=1,npoint(1)
         xx(i,1)=ampx(i)
c     xx(i,1)=(ampx(i)**2/re-re)/1.720787657571**2
c     yy(i,1)=-1.720787657571**2*damp(i)/(2*amp(i))
         yy(i,1)=-damp(i)/amp(i)
         write(51,*)xx(i,1),yy(i,1)
      end do
c     start=(start**2/re-re)/1.720787657571**2
c     end=(end**2/re-re)/1.720787657571**2
      heada=' Growth rate'
      call rita1a(xx,yy,start,end,ys,ye,npoint,1,npoint(1),'x',
     &     'alphai',heada,headb,headc,mbox,1,ihard)
      write(*,*)'again?'
      read(*,*)ans
      if (ans.eq.'y') go to 1

      end subroutine foustr
