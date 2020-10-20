c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine write_dx(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xl,zl,
     &     t,fileout,scalar)
c
c     Creates OpenDX readable files
c
      integer i,j,k,ll,nx,nyp,nz
      integer is,ie,js,je,scalar
      real xl,zl,t
      real ur(nx/2+1,nyp,nz,3+scalar),ui(nx/2+1,nyp,nz,3+scalar)
      real ar(nx/2+1,nyp,nz,2),ai(nx/2+1,nyp,nz,2)
      real dur(nx/2+1,nyp,nz,3,3),dui(nx/2+1,nyp,nz,3,3)
      real eta(nyp)
      character*80 fileout
      character*1000 string
      character*2 theta
      integer nfields
c
c     Complete domain
c
      is = 1
      ie = nx/2

      js = 1
      je = nyp
c
c     Some restrictions in output if desired
c     Note that j=1 is the top boundary, j=nyp the lower wall
c      is = 3500
c      ie = 4500

c      is = 1
c      ie = 1700/2+1
      
c      is = 1700/2
c      ie = 3400/2

c      js = 200
c      je = nyp

c
      write(*,*) 'x: ',is,' - ',ie,'(1 -',nx/2,')'
      write(*,*) 'y: ',js,' - ',je,'(1 -',nyp,')'
c
c     Write header file
c
      open(unit=22,file=trim(fileout)//'.general',form='formatted')
      write(22,'(a,a,a)') 'file = ',trim(fileout),'.dx'
      write(22,'(a,i4,a,i4,a,i4)') 'grid = ',
     &     (ie-is+1)*2,'x',je-js+1,'x',nz
      write(22,'(a)') 'format = msb ieee'
      write(22,'(a)') 'interleaving = record'
      write(22,'(a)') 'majority = column'
      write(22,'(a)') 'header = bytes 4'


      nfields = 3+scalar+3
c      nfields = 1
c
c     put nfields to 1 if only one field is written
c

      string = 'field = u'
      if (nfields.ge.2) then
         string = trim(string) // ',v'
      end if
      if (nfields.ge.3) then
         string = trim(string) // ',w'
      end if
      do i=1,scalar
         write(theta,'(a1,i1)') 't',i
         string = trim(string) // ',' // theta 
      end do
      if (nfields.gt.1) then
         string = trim(string) //',umean,lambda2,dudy'
      end if
      write(22,'(a,a)') trim(string)

      string = 'structure = scalar'
      do i=1,nfields-1
         string = trim(string) // ',scalar'
      end do
      write(22,'(a,a)') trim(string)
      
      string = 'type = float'
      do i=1,nfields-1
         string = trim(string) // ',float'
      end do
      write(22,'(a,a)') trim(string)
      
      string = 'dependency = positions'
      do i=1,nfields-1
         string = trim(string) // ',positions'
      end do
      write(22,'(a,a)') trim(string)


      write(22,'(a)') 'recordseparator =  bytes 8'
      write(22,'(a,f18.9,a)') 
     &     'positions = regular, irregular, regular, ',
     &  (is-1)*2*xl/nx ,','
      write(22,'(f18.9,a)') xl/nx,','
      do k=js,je
         if (k.eq.je) then
            write(22,'(f18.9)') eta(k)
         else
            write(22,'(f18.9,a)') eta(k),','
         end if
      end do
      write(22,'((f18.9,a,f18.9))') -zl/2.,',',zl /nz

      write(22,'(a)')
      write(22,'(a)') 'end'
      close(22)
c
c     Write meta file
c
      open(unit=22,file=trim(fileout)//'.meta',form='formatted')
      write(22,'(a20)') 'object 1 class field'
      write(22,'(a24,f18.9)') 'attribute "time" number ',t
      write(22,'(a3)') 'end'
      write(22,'(a2)') '  '
      close(22)
c
c     Write data file
c
c     The data structures is:
c     ur,ui(:,:,:,1:3) : velocities
c     ur,ui(:,:,:,4:3+scalar)  : scalars
c     ar,ai(:,:,:,1)   : udist
c     ar,ai(:,:,:,2)   : lambda2
c     dur,dui...
      open(unit=22,file=trim(fileout)//'.dx',form='unformatted')

      if (nfields.eq.1) then
         do ll=2,2
            write(22)
     &           ((( real(ar(i,j,k,ll),4),real(ai(i,j,k,ll),4),
     &           i=is,ie),j=js,je),k=1,nz)
         end do
         
      else

         do ll=1,3+scalar
            write(22)
     &           ((( real(ur(i,j,k,ll),4),real(ui(i,j,k,ll),4),
     &           i=is,ie),j=js,je),k=1,nz)
         end do
         
         do ll=1,2
            write(22)
     &           ((( real(ar(i,j,k,ll),4),real(ai(i,j,k,ll),4),
     &           i=is,ie),j=js,je),k=1,nz)
         end do
         
         write(22)
     &        ((( real(dur(i,j,k,1,2),4),real(dui(i,j,k,1,2),4),
     &        i=is,ie),j=js,je),k=1,nz)
      end if
      
      close(22)


      if (1.eq.0) then

c     PLOT3D format

      open(unit=22,file=trim(fileout)//'.xyz',form='unformatted',
     &     convert='big_endian')

c      write(22) 1
      write(22) (ie-is+1)*2,je-js+1,nz

      write(22) 
     &     (((real(xl/nx*(2*i-2),4),
     &        real(xl/nx*(2*i-1),4),
     &     i=is,ie),j=js,je),k=1,nz),
     &     (((real(eta(j),4),real(eta(j),4),
     &     i=is,ie),j=js,je),k=1,nz),
     &     (((real(-zl/2+zl/nz*(k-1),4),real(-zl/2+zl/nz*(k-1),4),
     &     i=is,ie),j=js,je),k=1,nz)

      close(22)


      open(unit=22,file=trim(fileout)//'.q',form='unformatted',
     &     convert='big_endian')

c      write(22,*) 1
      write(22) (ie-is+1)*2,je-js+1,nz
      write(22) real(1.,4),real(0.,4),real(100.,4),real(99.,4)

      do i=1,(nx/2+1)*nyp*nz
         dur(i,1,1,1,1) = ar(i,1,1,1)
         dui(i,1,1,1,1) = ai(i,1,1,1)

         dur(i,1,1,2,1) = ur(i,1,1,1)
         dui(i,1,1,2,1) = ui(i,1,1,1)

         dur(i,1,1,3,1) = ur(i,1,1,2)
         dui(i,1,1,3,1) = ui(i,1,1,2)

         dur(i,1,1,4,1) = ur(i,1,1,3)
         dui(i,1,1,4,1) = ui(i,1,1,3)

         dur(i,1,1,5,1) = ar(i,1,1,2)
         dui(i,1,1,5,1) = ai(i,1,1,2)
      end do

      write(22) 
     &     ((((real(dur(i,j,k,ll,1),4),real(dui(i,j,k,ll,1),4)
     &     ,i=is,ie),j=js,je),k=1,nz),ll=1,5)

     
      close(22)

      end if

      end subroutine write_dx
