c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine bfvort3(ur,ui,alfa,beta,my_node,prey)
c
c     Compute vorticities for 3d base flow
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,6),ui(memnx,memny,memnz,6)
      real alfa(nx/2*mbz),beta(nz)
      real app1(nyp),app2(nyp)
      real prey(nyp*2+15),w(nxp/2+1,nyp)
      integer my_node,i,x,z,nn,zbp,y
c
c     y-derivatives
c
      nn=0
      do i=1,3,2
         nn=nn+1
         do x=1,memnx
            do z=1,memnz
               do y=1,nyp
                  app1(y)=ur(x,y,z,i)
                  app2(y)=ui(x,y,z,i)
               end do
               call vchbf(app1,w,nyp,1,1,1,prey)
               call vchbf(app2,w,nyp,1,1,1,prey)
               call rdcheb(app1,nyp,1,1)
               call rdcheb(app2,nyp,1,1)
               call vchbb(app1,w,nyp,1,1,1,prey)
               call vchbb(app2,w,nyp,1,1,1,prey)
               do y=1,nyp
                  ur(x,y,z,7-i)=(-1)**nn*app1(y)*(2./real(nyp-1))
                  ui(x,y,z,7-i)=(-1)**nn*app2(y)*(2./real(nyp-1))
               end do
            end do
         end do
      end do
c
c     Vorticities
c
      do zbp=1,memnz
         z=my_node*memnz+zbp
         do x=1,memnx
            do y=1,nyp
c     omega_y
               ur(x,y,zbp,5)=-beta(z)*ui(x,y,zbp,1)+
     &              alfa(x)*ui(x,y,zbp,3)
               ui(x,y,zbp,5)= beta(z)*ur(x,y,zbp,1)-
     &              alfa(x)*ur(x,y,zbp,3)
c     omega_x
               ur(x,y,zbp,4)=ur(x,y,zbp,4)+beta(z)*ui(x,y,zbp,2)
               ui(x,y,zbp,4)=ui(x,y,zbp,4)-beta(z)*ur(x,y,zbp,2)
c     omega_z
               ur(x,y,zbp,6)=ur(x,y,zbp,6)-alfa(x)*ui(x,y,zbp,2)
               ui(x,y,zbp,6)=ui(x,y,zbp,6)+alfa(x)*ur(x,y,zbp,2)
            end do
         end do
      end do

      end subroutine bfvort3
