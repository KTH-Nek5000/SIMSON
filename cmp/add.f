c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine add(n,c,ur,ui,eta,fltype)
c
c     Adds field n to field 1 with coefficient c
c     If n=1 simply multiply by c
c
      implicit none

      include 'par.f'

      integer n,fltype
      real c
      real ur(memnx,memny,memnz,3,2),ui(memnx,memny,memnz,3,2)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real box2r(nx/2,mbz,nyp),box2i(nx/2,mbz,nyp)
      real eta(nyp)

      integer x,y,z,i,zb
      real um

      if (c.eq.1..and.n.eq.1) return
      do i=1,3
         do zb=1,nzc,mbz
            call getxy(boxr,boxi,zb,i,ur,ui)
            if (n.ne.1) then
               call getxy(box2r,box2i,zb,i+3,ur,ui)
            end if
            do z=zb,zb+mbz-1
               do y=1,nyp
                  if (z.eq.1.and.i.eq.1) then
                     um=eta(y)
                     if (fltype.eq.1.or.fltype.eq.4) then
                        um=box2r(1,1,1)+1.-eta(y)**2
                     else if (fltype.eq.3.or.fltype.eq.6.or.
     &                       fltype.eq.7.or.fltype.eq.8.or.
     &                       fltype.eq.9.or.fltype.eq.-1.or.
     &                       fltype.eq.-2) then
                        um=0.0
                     else if (fltype.eq.2.or.fltype.ge.5) then
                        um=box2r(1,1,1)+eta(y)
                        um = eta(y)
c                        um = 0.
                     else if (fltype.eq.-3) then
                        um = 0.
                     else
                        write(*,*) 'Flow type not supported by cmp',
     &                       fltype
                        stop
                     end if
c
c     If you don't want to subtract the parabolic base flow in
c     channel flow set um = 0.
c
                     if (n.ne.1) then
                        boxr(1,1,y)=boxr(1,1,y)+(box2r(1,1,y)-um)*c
                        boxi(1,1,y)=boxi(1,1,y)+box2i(1,1,y)*c
                        do x=2,nx/2
                           boxr(x,1,y)=boxr(x,1,y)+box2r(x,1,y)*c
                           boxi(x,1,y)=boxi(x,1,y)+box2i(x,1,y)*c
                        end do
                     else
                        boxr(1,1,y)=boxr(1,1,y)+(boxr(1,1,y)-um)*(c-1)
                        boxi(1,1,y)=boxi(1,1,y)*c
                        do x=2,nx/2
                           boxr(x,1,y)=boxr(x,1,y)*c
                           boxi(x,1,y)=boxi(x,1,y)*c
                        end do
                     end if
                  else
                     if (n.ne.1) then
                        do x=1,nx/2
                           boxr(x,z-zb+1,y) = boxr(x,z-zb+1,y)
     &                                       +box2r(x,z-zb+1,y)*c
                           boxi(x,z-zb+1,y) = boxi(x,z-zb+1,y)
     &                                       +box2i(x,z-zb+1,y)*c
                        end do
                     else
                        do x=1,nx/2
                           boxr(x,z-zb+1,y)=boxr(x,z-zb+1,y)*c
                           boxi(x,z-zb+1,y)=boxi(x,z-zb+1,y)*c
                        end do
                     end if
                  end if
               end do
            end do
            call putxy(boxr,boxi,zb,i,ur,ui)
         end do
      end do

      return

      end subroutine add
