c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wpgmr(pixm,nx,ny,mpix)
c
      integer nx,ny,mpix
      character*1 pixm(nx*ny+100)
c
      character*80 name
      character*80 line
      integer ifil,i
c
c     Finding a free file
c
      ifil=0
 1000 write(name,9100) ifil/10,mod(ifil,10)
 9100 format('plot',2i1,'.pgmr')
      open(unit=47,file=name,status='old',err=1010)
      close(unit=47)
      ifil=ifil+1
      goto 1000
 1010 close(unit=47)
      open(unit=47,file=name,status='new',access='direct',
     &recl=(nx*ny+100+7)/8*8,form='unformatted')
c
c     Write header
c
      pixm(1)='P'
      pixm(2)='5'
      pixm(3)=char(10)
      write(line,9000) nx,ny
 9000 format(i14,'  ',i14)
      do i=1,30
         pixm(3+i)=line(i:i)
      end do
      pixm(34)=char(10)
      write(line,9010) mpix
 9010 format(50x,i15)
      do i=1,65
         pixm(34+i)=line(i:i)
      end do
      pixm(100)=char(10)
      write(47,rec=1) pixm
      close(unit=47)

      return

      end subroutine wpgmr
