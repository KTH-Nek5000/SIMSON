      subroutine rwavesbl(omega,ucwa,uswa,vcwa,vswa,ndxwa,
     &           dstar,my_node)

      implicit none

      include 'par.f'

      real ucwa(nyp,180),uswa(nyp,180),vcwa(nyp,180),vswa(nyp,180)
      real omega,dstar,ucmax,usmax,vcmax,vsmax,fmax
      integer i,n,ndxwa,my_node

      ucmax=0.
      usmax=0.
      vcmax=0.
      vsmax=0.
c
c     Reads wave from the file waves.d
c
      open(unit=10,status='old',file='waves.d')
      read(10,*) omega,ndxwa
c
c     Rescale in channel coordinates
c
      omega=omega/dstar
      if (ndxwa.gt.180) then
         write(*,*) 'nxdwa too large for ucwa etc. in rwavesbl.f'
         stop
      end if
      do i=1,ndxwa
         do n=1,nyp
            read(10,*) ucwa(n,i),uswa(n,i)
            if (abs(ucwa(n,i)).gt.ucmax) ucmax=abs(ucwa(n,i))
            if (abs(uswa(n,i)).gt.usmax) usmax=abs(uswa(n,i))
         end do
      end do

      do i=1,ndxwa
         do n=1,nyp
            read(10,*) vcwa(n,i),vswa(n,i)
            if (abs(vcwa(n,i)).gt.vcmax) vcmax=abs(vcwa(n,i))
            if (abs(vswa(n,i)).gt.vsmax) vsmax=abs(vswa(n,i))
         end do
      end do
      close(unit=10)

      if (ucmax.gt.fmax) fmax=ucmax
      if (usmax.gt.fmax) fmax=usmax
      if (vcmax.gt.fmax) fmax=vcmax
      if (vsmax.gt.fmax) fmax=vsmax

      do n=1,nyp
         do i=1,ndxwa
            ucwa(n,i)=ucwa(n,i)/fmax
            uswa(n,i)=uswa(n,i)/fmax
            vcwa(n,i)=vcwa(n,i)/fmax
            vswa(n,i)=vswa(n,i)/fmax
         end do
      end do

      if (my_node.eq.0) then
         write(*,*) 'Reading waves with omega:',omega*dstar
      end if

      end subroutine rwavesbl
