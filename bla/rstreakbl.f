c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rstreakbl(betast,omegast,ndxst,
     &     uust_r,uust_i,vvst_r,vvst_i,wwst_r,wwst_i,str_nam1,str_nam2)
c
c     Read streak data from files str_nam1 and str_nam2
c
      implicit none

      include 'par.f'

      real uust_r(nyp,180,2),vvst_r(nyp,180,2),wwst_r(nyp,180,2)
      real uust_i(nyp,180,2),vvst_i(nyp,180,2),wwst_i(nyp,180,2)
      real betast(2),omegast(2)
      integer ndxst,ndxst1
      character*80 str_nam1,str_nam2
      real phase_r,phase_i,app

      integer i,n

      open(unit=10,status='old',file=str_nam1)

      read(10,*) betast(1),omegast(1),ndxst
c
c     Note that betast,omegast are already in code scaling
c
      if (ndxst.gt.180) then
         write(ioe,*) 'not enough streamwise positions in rstreakbl'
         call stopnow(4334543)
      end if

      do i=1,ndxst
         do n=1,nyp
            read(10,*) uust_r(n,i,1),uust_i(n,i,1),
     &           vvst_r(n,i,1),vvst_i(n,i,1),
     &           wwst_r(n,i,1),wwst_i(n,i,1)

         end do
      end do
      read(10,*) phase_r,phase_i
      close(unit=10)

      do i=1,ndxst
         do n=1,nyp
            app=uust_r(n,i,1)*phase_r-uust_i(n,i,1)*phase_i
            uust_i(n,i,1)=uust_i(n,i,1)*phase_r+uust_r(n,i,1)*phase_i
            uust_r(n,i,1)=app
            app=vvst_r(n,i,1)*phase_r-vvst_i(n,i,1)*phase_i
            vvst_i(n,i,1)=vvst_i(n,i,1)*phase_r+vvst_r(n,i,1)*phase_i
            vvst_r(n,i,1)=app
            app=wwst_r(n,i,1)*phase_r-wwst_i(n,i,1)*phase_i
            wwst_i(n,i,1)=wwst_i(n,i,1)*phase_r+wwst_r(n,i,1)*phase_i
            wwst_r(n,i,1)=app
         end do
      end do

      open(unit=10,status='old',file=str_nam2)

      read(10,*) betast(2),omegast(2),ndxst1

      if (ndxst.ne.ndxst1) then
         write(ioe,*) 'number of streak positions differ'
         call stopnow(5543543)
      end if

      do i=1,ndxst
         do n=1,nyp
            read(10,*) uust_r(n,i,2),uust_i(n,i,2),
     &           vvst_r(n,i,2),vvst_i(n,i,2),
     &           wwst_r(n,i,2),wwst_i(n,i,2)
         end do
      end do
      read(10,*)phase_r,phase_i
      close(unit=10)

      do i=1,ndxst
         do n=1,nyp
            app=uust_r(n,i,2)*phase_r-uust_i(n,i,2)*phase_i
            uust_i(n,i,2)=uust_i(n,i,2)*phase_r+uust_r(n,i,2)*phase_i
            uust_r(n,i,2)=app
            app=vvst_r(n,i,2)*phase_r-vvst_i(n,i,2)*phase_i
            vvst_i(n,i,2)=vvst_i(n,i,2)*phase_r+vvst_r(n,i,2)*phase_i
            vvst_r(n,i,2)=app
            app=wwst_r(n,i,2)*phase_r-wwst_i(n,i,2)*phase_i
            wwst_i(n,i,2)=wwst_i(n,i,2)*phase_r+wwst_r(n,i,2)*phase_i
            wwst_r(n,i,2)=app
         end do
      end do

      end subroutine rstreakbl
