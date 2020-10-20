c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program hdf2bla

      implicit none
c
c     Read in HDF velocity field and writes 'bla' .u fltype=6 file
c
#ifdef HDF
      include 'hdf.inc'
#endif

      character*20 filein, fileout, fbase
      integer i,j,k,ll

      integer rank, mode
      integer status, sd_id, sfstart, sds_id, sfrattr,sfrdata,sfendacc
      integer sfselect,sfginfo,dsnum,nsets
      integer attr, type, sfend
      character*80 sds_name

      real re,xl,yl,zl,xs,dstar,t,pi
      real bstart, blength
      integer fltype, nx2,nyp2,nz2,nx,nyp,nz,nfzsym
      logical pou
      real zs,argx,argz,hr

      integer,allocatable :: stride(:), start(:),dimsizes(:)

      real,allocatable :: data3D(:,:,:,:),data2D(:,:,:)
      real,allocatable :: ur(:,:,:,:),ui(:,:,:,:),urx(:)
      real,allocatable :: alfa(:),beta(:)
      real,allocatable :: prexn(:),prezn(:),wr(:),wi(:)
      real,allocatable :: cx(:),sx(:),sa(:),ca(:)

      pi=4.*atan(1.)

#ifndef HDF
      write(*,*) 'HDF tools support not activated. Recompile with',
     &           ' proper options.'
#else
      write(*,*) 'This program converts velocity fields stored on',
     &           ' HDF-format to bla readable format (fltype=6).'
      write(*,*)
      write(*,*) 'usage: hdf2bla  [inputfile baseflow outputfile]'

      call getarg(1,filein)
      call getarg(2,fbase)
      call getarg(3,fileout)

      write(*,*)
      write(*,*) 'Read from file : ',trim(filein)
      write(*,*) 'Read baseflow from file : ',trim(fbase)
      write(*,*) 'Write to file : ',trim(fileout)
c
c     Open HDF-file
c
      sd_id = sfstart(filein,DFACC_READ)
      if (sd_id .eq. -1) then
         write(*,*) 'HDF file does no exists'
         goto 9000
      end if

      nsets = dsnum(filein)
      write(*,*)

      sds_id = sfselect(sd_id,1)
      write(*,*) 'Give rank (dimensions): '
      read(*,*) rank

      allocate(dimsizes(rank),stride(rank), start(rank))
c
c     Get info of the HDF-file
c
      status = sfginfo(sds_id, sds_name, rank, dimsizes, type, attr)

      if (rank.eq.2) then
         nx = dimsizes(1)
         nyp = dimsizes(2)
         nz = 4
         write(*,*) 'nz is chosen to: ', nz
         allocate(data2D(nx,nyp,rank))
      else if (rank.eq.3) then
         nx = dimsizes(1)
         nyp = dimsizes(2)
         nz = dimsizes(3)
         allocate(data3D(nx,nyp,nz,rank))
      else
         write(*,*) 'Only Rank 2 and 3 (2D and 3D)'
         goto 9000
      end if

      do i=1,rank
         stride(i) = 1
         start(i) = 0
      end do

      write(*,*) 'Give the dataset to start to read:'
      read(*,*)  mode

      if (rank.eq.2) then
         j = 2*mode-3
         do i=1,2
            sds_id = sfselect(sd_id,j+i)
            status = sfrdata(sds_id, start, stride, dimsizes
     &           , data2D(1,1,i))
            if (status.eq.-1) then
               write(*,*)  'Reading failed'
               goto 9000
            end if
            status = sfendacc(sds_id)
         end do
         status = sfend(sd_id)
      else if (rank.eq.3) then
         j = 3*mode-4
         do i=1,3
            sds_id = sfselect(sd_id,j+i)
            status = sfrdata(sds_id, start, stride, dimsizes
     &           ,data3D(1,1,1,i))
            if (status.eq.-1) then
               write(*,*)  'Reading failed'
               goto 9000
            end if
            status = sfendacc(sds_id)
         end do
         status = sfend(sd_id)
      end if
      deallocate(stride,start,dimsizes)
c
c     Write data into bla-format
c
      allocate(ur(nx/2+1,nyp,nz,3),ui(nx/2+1,nyp,nz,3))
      write(*,*) 'Writing to bla-format'

      do ll=1,3
         do k=1,nz
            do j=1,nyp
               do i=1,nx/2
                  if (rank.eq.2) then
                     ur(i,j,k,ll) = data2D(2*i-1,j,ll)
                     ui(i,j,k,ll) = data2D(2*i,j,ll)
                     if (ll.eq.3) then
                        ur(i,j,k,ll) = 0.
                        ui(i,j,k,ll) = 0.
                     end if
                  else if (rank.eq.3) then
                     ur(i,j,k,ll) = data3D(2*i-1,j,k,ll)
                     ui(i,j,k,ll) = data3D(2*i,j,k,ll)
                  end if
               end do
               ur(nx/2+1,j,k,ll) = 0.
               ui(nx/2+1,j,k,ll) = 0.
               if (k.eq.nz/2+1) then
                  do i=1,nx/2
                     ur(i,j,k,ll) = 0.
                     ui(i,j,k,ll) = 0.
                  end do
               end if
            end do
         end do
      end do

c
c     Read in the base flow from bla.file
c
      open(unit=12,file=fbase,status='old',form='unformatted')
      read(12) re,pou,xl,zl,t,xs
      read(12) nx2,nyp2,nz2
      read(12) fltype,dstar
      read(12) bstart,blength
      close(12)

      write(*,*) 'Reading base-flow'
      re = re * dstar
      xl = xl / dstar
      zl = zl / dstar
      yl = 2. / dstar
      t  = t  / dstar

      write(*,*)
      write(*,*) 'Re             : ',re
      write(*,*) 't              : ',t
      write(*,*) 'Resolution     : ',nx,nyp,nz
      write(*,'(A,3F12.5)') ' Box            : ',xl,yl,zl
c
c     Compute alfa and beta
c
      allocate(alfa(nx/2),beta(nz))
      do i=1,nx/2
         alfa(i) = 2.*pi/xl*real(i-1)
      end do
      beta(1) = 0.
      do k=2,nz/2+1
         beta(k) = 2.*pi/zl*real(k-1)
         beta(nz+2-k) = -2.*pi/zl*real(k-1)
      end do
c
c     Do forward FFT to Fourier space
c
      write(*,*) 'Do FFT...'
      allocate(prexn(nx+25))
      allocate(prezn(nz*2+25))
      allocate(wr((nx/2+1)*nyp*nz),wi((nx/2+1)*nyp*nz))
      call vrffti(nx,prexn,i)
      call vcffti(nz,prezn,i)

      do ll=1,3
         do j=1,nyp
            call vcfftf(ur(1,j,1,ll),ui(1,j,1,ll),wr,wi,
     &           nz,nx/2,(nx/2+1)*nyp,1,prezn)
            call vrfftf(ur(1,j,1,ll),ui(1,j,1,ll),wr,wi,
     &           nx,nz,1,(nx/2+1)*nyp,prexn)
         end do
      end do
c
c     Shift back data
c
      write(*,*) 'Shift back data...'
      xs = -xl/2.
      zs = 0.
      allocate(cx(nx/2),sx(nx/2),ca(nx/2),sa(nx/2))
      do i=1,nx/2
         argx = -xs*alfa(i)
         cx(i) = cos(argx)
         sx(i) = sin(argx)
      end do
      do k=1,nz
         argz = -zs*beta(k)
         do i=1,nx/2
            ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz)
            sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz)
         end do
         do ll=1,3
            do j=1,nyp
               do i=1,nx/2
                  hr=ur(i,j,k,ll)*ca(i)-ui(i,j,k,ll)*sa(i)
                  ui(i,j,k,ll)=ui(i,j,k,ll)*ca(i)+ur(i,j,k,ll)*sa(i)
                  ur(i,j,k,ll)=hr
               end do
            end do
         end do
      end do
      deallocate(cx,sx,ca,sa)
c
c     Write data to disc in bla-format and normalize (fft)
c
      allocate(urx(nx))
      write(*,*) 'write to disc: ', fileout
      re = re / dstar
      xl = xl * dstar
      zl = zl * dstar
      yl = 2. * dstar
      t  = t  * dstar

      open(unit=11,file=fileout,form='unformatted')
      write(11) re,pou,xl,zl,t,0.
      write(11) nx,nyp,nz,nfzsym
      write(11) fltype,dstar
      write(11) bstart,blength,0.,0.

      do ll=1,3
         do k=1,nz
            do j=1,nyp
               do i=1,nx/2
                  urx(2*i-1) = ur(i,j,k,ll)/real(nx*nz)
                  urx(2*i) = ui(i,j,k,ll)/real(nx*nz)
               end do
               write(11) (urx(i),i=1,nx)
            end do
         end do
      end do
      deallocate(urx)
      deallocate(ur,ui)

      call flush(11)
      close(unit=11)
 9000 continue
#endif

c
c     Finalize subroutine
c
      end program hdf2bla
