c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine foubranch(plane,eta,xr,xl,heada,headb,re,mbox,ihard)
c
c     Calculates branch I,II for the neutral stability
c     curve and plots the curve of the maximum over y
c     of the choosen variabl.
c
c     Originally by Stellan Berlin
c
      implicit none

      include 'par.f'

      character*1 ans
      integer mbox,ihard,deg,ans2
      integer i,j,x,y,nex,nbr,npoint(1),xsp
      real plane(0:nx,nyp),eta(nyp),xr,xl,re,amp(4*nx),damp(4*nx)
      real ampx(4*nx),coeff(10)
      real u(nyp),y2(nx+nyp),xex(nx+nyp),yex(nx+nyp),umax(nx),ymax(nx)
      real ux(nx),xbr(nx),ybr(nx),lnu(4*nx),xx(4*nx,1),yy(4*nx,1)
      real start,end,norm,ys,ye
      character*80 heada,headb,headc
c
c     Find umax(y)
c
      do x=1,nx
         do y=1,nyp
            u(y)=plane(x,y)
c            if (x.ge.130.and.x.le.150) write(50,*) x,y,u(y)
         end do
         do y=1,nyp
            yex(y)=0.
            xex(y)=0.
         end do
         call spline(eta,u,y2,nyp,xex,yex,nex)
         ymax(x)=0.
         umax(x)=0.
         do i=1,nex
            if (abs(umax(x)).lt.abs(yex(i))) then
               umax(x)=yex(i)
               ymax(x)=xex(i)
            end if
         end do
c         write(51,*) x,umax(x),ymax(x)
      end do

      do x=1,nx
         ux(x)=real(x-nx/2-1)/real(nx)*xl+xr
         ux(x)=1.720787657571*sqrt(re*(ux(x)+re/1.720787657571**2))
c         xx(x,1)=ux(x)
c         yy(x,1)=umax(x)
      end do
c      npoint(1)=nx
c      write(headc,9601) xbr(1),xbr(2)
c 9601 format('Branch I at Re1=',F6.2,' Branch II at Re1=',F6.2)

c      call rita1a(xx,yy,1.,1.,1.,1.,npoint,1,nx,'Re1',
c     &    'logA/Ao',heada,headb,headc,1,1,-1)
c      read(*,9701) ans
c 9701 format(a1)

      write(*,*) 'Give first real plot point'
      read(*,*) xsp
      do x=1,xsp-1
         umax(x)=umax(xsp)+(xsp-x)*(umax(xsp)-umax(xsp+1))/2
      end do
c
c     Make amplitude curve and find branch I,II
c
      call spline(ux,umax,y2,nx,xbr,ybr,nbr)

      do x=1,nx
         xx(x,1)=0
         yy(x,1)=0
      end do
      do x=1,nx
         xx(x,1)=ux(x)
         yy(x,1)=umax(x)
      end do

      npoint(1)=nx
      write(headc,'(a)') '    '
      call rita1a(xx,yy,1.,1.,1.,1.,npoint,1,nx,'Re1',
     &     'Amplitunde',heada,headb,headc,mbox,1,ihard)

      if (ybr(1).eq.0.or.ybr(2).eq.0) then
         ans='y'
      else
         write(*,*) 'Change default options ? (y/n)'
         read(*,9100) ans
 9100    format(a1)
      end if

      npoint(1)=nx

      if (ans.eq.'Y'.or.ans.eq.'y') then
         write(*,9200) ux(1),ux(nx)
 9200    format('Give start and end for re1 (',f4.0,'-',f4.0,')')
         read(*,*) start,end
         if (start.lt.ux(1).or.start.gt.ux(nx)) start=ux(1)
         if (end.gt.ux(nx).or.end.lt.ux(1)) end=ux(nx)
         write(*,*) 'Give ystart and yend '
         read(*,*) ys,ye
         write(*,9300) nx,4*nx
 9300    format('Give number of re1 points (0 for default ',i4,
     &        ' max ',i4,')')
         read(*,*) npoint(1)
         if (npoint(1).eq.0) npoint(1)=nx
 5       write(*,9350) ybr(1)
 9350    format('Give normalisation factor (0 for default ',e8.3,')')
         read(*,*) norm
         if (norm.eq.0) norm=ybr(1)
         if (norm.le.0) then
            write(*,*) 'The normalisation factor must be > 0'
            goto 5
         end if
         write(*,*) ' Do you want curve fitting ? '
         write(*,*) ' 0  No curve fitting'
         write(*,*) '>0  Degree of polynomial'
         read(*,*) deg
      else
         start=ux(1)
         do i=nx,1,-1
            if (ux(i).lt.xbr(2)+70) goto 10
         end do
 10      end=ux(i)
         norm=ybr(1)
         ys=1.
         ye=1.
      end if

      do i=1,npoint(1)
         ampx(i)=start+i*(end-start)/npoint(1)
      end do

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
      else
         call splint(ux,umax,y2,nx,ampx,amp,damp,npoint(1))
      end if

      if (amp(npoint(1)).lt.0.) then
         write(*,*) '    '
         write(*,*) '**  Try with lower degre of polynomial  **'
         write(*,*) '    '
         write(*,*) '    '
         write(*,*) '    '
         go to 22
      end if
      do i=1,npoint(1)
         lnu(i)=log(amp(i)/norm)
c         write(51,*) ampx(i),amp(i)
      end do

      do i=1,npoint(1)
         xx(i,1)=0
         yy(i,1)=0
      end do

      do i=1,npoint(1)
         xx(i,1)=ampx(i)
         yy(i,1)=lnu(i)
      end do

      write(headc,9400) xbr(1),xbr(2)
 9400 format('Branch I at Re1=',F6.2,' Branch II at Re1=',F6.2)

      call rita1a(xx,yy,start,end,ys,ye,npoint,1,npoint(1),'Re1',
     &    'logA/Ao',heada,headb,headc,mbox,1,ihard)

      write(*,*) 'Do you want to see the growth rate curve (y/n)'
      read(*,9500) ans
 9500 format(a1)

      if (ans.eq.'Y'.or.ans.eq.'y') then
         do i=1,npoint(1)
            xx(i,1)=0
            yy(i,1)=0
         end do

         do i=1,npoint(1)
            xx(i,1)=ampx(i)
            yy(i,1)=-1.720787657571**2*damp(i)/(2*amp(i))
         end do

         call rita1a(xx,yy,start,end,ys,ye,npoint,1,npoint(1),'Re1',
     &       'alphai',heada,headb,headc,mbox,1,ihard)

         if (deg.gt.0) then
            write(*,*) 'Do you want to calculate branch I,II ?'
            write(*,*) '0  no'
            write(*,*) '1  yes'
            write(*,*) '2  yes with my initial guesses'
            read(*,*) ans2
            if (ans2.gt.0) then
               if (ans2.gt.1) then
                  write(*,*) 'give initial guesses'
                  read(*,*) xbr(1),xbr(2)
               end if
               call newton(xbr(1),coeff,deg)
               call newton(xbr(2),coeff,deg)
               write(headc,9400) xbr(1),xbr(2)
               call rita1a(xx,yy,start,end,ys,ye,npoint,1,npoint(1),
     &              'Re1','alphai',heada,headb,headc,mbox,1,ihard)
            end if
         end if
      end if

 22   return

      end subroutine foubranch
