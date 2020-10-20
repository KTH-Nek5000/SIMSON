c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wraw(pixm,nx,ny,ifil)
c
c     Opens a raw data file and writes nx*ny bytes
c
      implicit none

      integer nx,ny,ifil
      character*1 pixm(nx*ny)

      character*80 name
      integer darecl
c
c     Finding a free file
c
 1000 write(name,9100) ifil/100,mod(ifil,100)/10,mod(ifil,10)
 9100 format('plot',3i1,'.raw')
      open(unit=47,file=name,status='old',err=1010)
      close(unit=47)
      ifil=ifil+1
      goto 1000
 1010 close(unit=47)
c      call frecl((nx*ny+7)/8*8,darecl)
      darecl = 1
      open(unit=47,file=name,status='new',access='direct',
     &     recl=darecl,form='unformatted')
      write(47,rec=1) pixm
      close(unit=47)
      return

      end subroutine wraw
