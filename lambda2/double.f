c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program double

      implicit none

      integer,parameter :: scalar=2

      real,allocatable :: urx(:)

      real re,xl,yl,zl,t,xs,rlam,dstar,spanv,bstart,blength,pr(scalar)
      real pi,umeanr,umeani,m1(scalar)
      real zs,argx,argz,hr

      integer i,j,k,ll,kk
      integer nx,nyp,nz,fltype

      logical pou,cb

      character(len=80) filein,fileout

      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '*                             Simson        '//
     &     '                      *'
      write(*,*) '*                        lambda2 $Rev$'//
     &     '                      *'
      write(*,*) '*                                           '//
     &     '                      *'
      write(*,*) '********************************************'//
     &     '***********************'
      write(*,*)

      write(*,*) 'usage: double  [inputfile outputfile]'


      call getarg(1,filein)
      call getarg(2,fileout)

      if (len(trim(filein)).eq.0) then
         filein = 'input.u'
      end if

      if (len(trim(fileout)).eq.0) then
         fileout = 'output.u'
      end if

      write(*,*)
      write(*,*) 'Read from file : ',trim(filein)
      write(*,*) 'Write to files : ',trim(fileout)

      open(unit=12,file=filein,status='old',form='unformatted')
      open(unit=13,file=fileout,status='unknown',form='unformatted')
c
c     Read meta data
c
      read(12) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)

      read(12) nx,nyp,nz
      fltype= 0
      dstar = 1.0
      read(12,err=1010) fltype,dstar
 1010 continue
c
c     Read additional info for specific flow types
c
      rlam = 0.
      spanv = 0.
      if (fltype.lt.0) read(12) rlam
      if (fltype.eq.6) read(12) bstart,blength
      if (fltype.ge.7) read(12) bstart,blength,rlam,spanv

      write(13) re,pou,xl,2*zl,t,xs,(pr(i),m1(i),i=1,scalar)
      write(13) nx,nyp,2*nz,0
      write(13) fltype,dstar
      if (fltype.lt.0) write(13) rlam
      if (fltype.eq.6) write(13) bstart,blength
      if (fltype.ge.7) write(13) bstart,blength,rlam,spanv

      allocate(urx(nx))
      do ll=1,3+scalar
         do k=1,nz

            do j=1,nyp
               read(12) (urx(i),i=1,nx)
               write(13) (urx(i),i=1,nx)
            end do
            urx = 0.
            do j=1,nyp
               write(13) (urx(i),i=1,nx)
            end do

         end do
      end do
      deallocate(urx)

      end program double
