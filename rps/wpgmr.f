c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wpgmr(pixm,nx,ny,ifil)
c
c     Opens a pgmr data file and writes nx*ny bytes
c
      implicit none

      integer i,nx,ny,ifil
      character*1 pixm(nx*ny)
      character*1 data(nx*ny+11)

      character*80 name
      character*15 string
      integer darecl
c
c     Finding a free file
c
 1000 write(name,9100) ifil/1000,
     &     mod(ifil,1000)/100,mod(ifil,100)/10,mod(ifil,10)
 9100 format('plot',4i1,'.pgm')
      open(unit=47,file=name,status='old',err=1010)
      close(unit=47)
      ifil=ifil+1
      goto 1000
 1010 close(unit=47)
      call frecl((nx*ny+15+7)/8*8,darecl)
      open(unit=47,file=name,status='new',access='direct',
     &     recl=darecl,form='unformatted')

      write(string, '(a3,i3,a1,i3,a1,i3,a1)')
     &     'P5 ',nx,' ',ny,' ',255,' '

      do i=1,15
         data(i) = string(i:i)
      end do
      do i=16,15+nx*ny
         data(i) = pixm(i-15)
      end do

      write(47,rec=1) data
      close(unit=47)

      return

      end subroutine wpgmr
