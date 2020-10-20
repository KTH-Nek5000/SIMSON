c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
c
c     setenv F_UFMTENDIAN "big;little:12"
c     unsetenv F_UFMTENDIAN !!!
      program lambda2
c
c     Data structure:
c     ur,ui (3+scalar)
c     dur,dui (9)
c     ar,ai (1)
c
c     For velocities: 13 velocity fields, i.e. about 4.3*disk space 
c
      implicit none
      integer,parameter :: scalar=0
c
c     Variables read from bla file with fixed 8-byte accuracy
c
      real*8,allocatable :: urx(:)
      real*8 re,xl,yl,zl,t,xs,rlam,dstar,spanv,bstart,blength,pr(scalar)
      real*8 m1(scalar),gr(scalar)
c
c     The following variables can be either 4 or 8 bytes
c
      integer,external :: omp_get_thread_num,omp_get_num_threads
      integer :: itot


      real,allocatable :: ur(:,:,:,:),ui(:,:,:,:)
      real,allocatable :: omr(:,:,:,:),omi(:,:,:,:)
      real,allocatable :: ar(:,:,:,:),ai(:,:,:,:)
      real,allocatable :: dur(:,:,:,:,:),dui(:,:,:,:,:)

      real,allocatable :: alfa(:),beta(:),eta(:)
      real,allocatable :: cx(:),sx(:),sa(:),ca(:)

      real,allocatable :: prexn(:),prezn(:),prey(:)
      real,allocatable :: wr(:,:),wi(:,:),w3(:,:),xcoord(:),zcoord(:)
      
      real,allocatable :: scalings(:,:)

      real pi,umeanr,umeani
      real zs,argx,argz,hr
c
c     Clenshaw stuff
c
      integer nypc
      real dy,ymax,ymin
      real,allocatable :: urc(:,:,:,:),uic(:,:,:,:)
      real,allocatable :: omrc(:,:,:,:),omic(:,:,:,:)
      real,allocatable :: arc(:,:,:,:),aic(:,:,:,:)

      character ch(4)

      integer i,j,k,ll,kk,ii,mm,jj
      integer nx,nyp,nz,fltype,ivort
      integer nxin,nypin,nzin
      real uminr,umaxr,umini,umaxi
c
c     Various options
c
      logical pou,cb,vort,uniform_y,subset
      integer xfac, zfac, numvar, numscalars
      integer, allocatable :: vars(:)
      integer, allocatable :: scalars(:)
      integer :: bounds(6)
      integer wu,wv,ww,wl2,wo1,wo2,wo3,wo4,wf,ws

      character(len=80) filein,fileout,inputfile
      character(len=2) type
      character dummy

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

      pi=4.*atan(1.)
#ifdef HDF
      write(*,*) 
     &     'usage: lambda2  [inputfile outputfile] [dx|ef|ec|vt|df]'
#else
      write(*,*) 
     &     'usage: lambda2  [inputfile outputfile settingsfile]'
#endif

      call getarg(1,filein)
      call getarg(2,fileout)
      call getarg(3,inputfile)

      if (len(trim(filein)).eq.0) then
         filein = 'input.u'
      end if

      if (len(trim(fileout)).eq.0) then
         fileout = 'output'
      end if

      write(*,*)
      write(*,*) 'Note that the number of scalar variables present'
      write(*,*) 'in the .u file needs to be correctly given in'
      write(*,*) 'the lambda2.f file.'
      write(*,*)
      write(*,*) 'Run lambda2 without environment variable'
      write(*,*) 'export F_UFMTENDIAN='
      write(*,*)
      write(*,*) 'Read from file : ',trim(filein)
#ifdef HDF
      write(*,*) 'Write to files : ',trim(fileout),
     &     '(.general|.dx|.meta) or (.case|.geo|.vel|.l2) or (.hdf)'
#else
      write(*,*) 'Write to files : ',trim(fileout),
     &     '(.general|.dx|.meta) or (.case|.geo|.vel|.l2)'
#endif

      write(*,*)
      write(*,*) 'Accuracy of Simson input file  : ',kind(urx)
      write(*,*) 'Internal accuracy              : ',kind(ur)

      if (len(trim(inputfile)).ne.0) then
         open(unit=10,status='old',file=trim(inputfile))
         read(10,9000) type
 9000    format(a2)
         write(*,*) 'Output format                  : ',trim(type)

         read(10,*) xfac
         write(*,*) 'Resolution factor in x (post)  : ',xfac

         read(10,*) zfac
         write(*,*) 'Resolution factor in z (post)  : ',zfac

         read(10,*) vort
         if (vort) then
            write(*,*) 'Computing vorticity'
         else
            write(*,*) 'Not computing vorticity'
         end if

         read(10,*) uniform_y
         if (uniform_y) then
            write(*,*) 'Uniform grid in y (Clenshaw interpolation)'
         end if

         if (uniform_y) then
            read(10,*) nypc
            read(10,*) ymin
            read(10,*) ymax

            write(*,*) '   Uniform grid points      : ',nypc
            write(*,*) '   Box y-range lower bound  : ',ymin
            write(*,*) '   Box y-range upper bound  : ',ymax
         else
c
c     Read dummy lines
c
            read(10,'(A1)') dummy
            read(10,'(A1)') dummy                   
            read(10,'(A1)') dummy              
         end if
c
c     Write out only a fraction of the box
c     (this can already be done for uniform grid in y)
c
         read(10,*) subset
         if (subset) then
            write(*,*) 'Writing a fraction of the box, indices:'
            read(10,*) bounds(1)
            read(10,*) bounds(2)
            read(10,*) bounds(3)
            read(10,*) bounds(4)
            read(10,*) bounds(5)
            read(10,*) bounds(6)
            write(*,*) '   x: ',bounds(1),' - ',bounds(2)
            write(*,*) '   y: ',bounds(3),' - ',bounds(4)
            write(*,*) '   z: ',bounds(5),' - ',bounds(6)
         else
c
c     Read dummy lines
c
            read(10,'(A1)') dummy
            read(10,'(A1)') dummy                   
            read(10,'(A1)') dummy              
            read(10,'(A1)') dummy         
            read(10,'(A1)') dummy              
            read(10,'(A1)') dummy
         end if
c
c     Figure out which variables are written out
c
         read(10,*) numvar  
         allocate(vars(numvar))
         do i=1,numvar
            read(10,*) vars(i)
            if (vars(i) .eq. 1) then
               write(*,*) 'Writing out u'
               wu=1
            elseif (vars(i) .eq. 2) then
               write(*,*) 'Writing out v'
               wv=1
            elseif (vars(i) .eq. 3) then
               write(*,*) 'Writing out w'
               ww=1
            elseif (vars(i) .eq. 4) then
               write(*,*) 'Writing out lambda2'
               wl2=1
            elseif (vars(i) .eq. 5) then
               write(*,*) 'Writing out omega_x'
               wo1=1
            elseif (vars(i) .eq. 6) then
               write(*,*) 'Writing out omega_y'
               wo2=1
            elseif (vars(i) .eq. 7) then
               write(*,*) 'Writing out omega_z'
               wo3=1
            elseif (vars(i) .eq. 8) then
               write(*,*) 'Writing out ufluct'
               wf=1
            elseif (vars(i) .eq. 9) then
               write(*,*) 'Writing out vorticity'
               wo4=1
            end if
         end do
c
c     Figure out which scalars (if any) are written out
c
         if (scalar.gt.0) then
            read(10,*) numscalars 
            allocate(scalars(numscalars))
            do i=1,numscalars
               read(10,*) scalars(i)      
               write(*,*) 'Writing out scalar number',scalars(i)
            end do
         end if

         close(10)
      else
         write(*,*) '****************************************'
         write(*,*) 'No input file, using default parameters.'
         write(*,*) '****************************************'
         type = 'vt'
         xfac = 1
         zfac = 1
         subset = .false.
         uniform_y = .false.
         vort = .false.
         ymax = yl
         numvar = 3
         numscalars = 0
         allocate(scalars(numscalars))
         allocate(vars(numvar))
         vars(1) = 1
         vars(2) = 4
         vars(3) = 8
      end if



c
c     Testing openmp
c
      itot = 1
!$omp parallel
      itot = omp_get_num_threads()
!$omp end parallel
      write(*,*) 'Running on ',itot,' threads.'

c
c     Probing little/big endian and opening file
c
      open(unit=12,file=filein,status='old',access='direct',recl=4,
     &     form='unformatted')
      read(12,rec=1) ch
      close(12)

      write(*,*) '===>',ch,'<===',
     &     ichar(ch(1)),ichar(ch(2)),ichar(ch(3)),ichar(ch(4))

      if (ichar(ch(1)).eq.0) then
c     big endian
         write(*,*) 'Open big endian file'
         open(unit=12,file=filein,status='old',form='unformatted',
     &        convert='big_endian')
      else
c     little endian
         write(*,*) 'Open little endian file'
         open(unit=12,file=filein,status='old',form='unformatted',
     &        convert='little_endian')
      end if
c         open(unit=12,file=filein,status='old',form='unformatted')

c
c     Read meta data
c
      rewind(12)
c      write(*,*) scalar
      read(12) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)
      read(12) nxin,nypin,nzin
      fltype= 0
      dstar = 1.0
      read(12,err=1010) fltype,dstar

      write(*,*) dstar

c
c     Description of parameters:
c     re      Reynolds number (rescaled, real: re*dstar)
c     pr      Prandtl number (only read for scalar=1)
c     pou     Not in use anymore, always set to false
c     xl      streamwise length (rescaled, real: xl/dstar)
c     zl      spanwise width (rescaled, real: xl/dstar)
c     t       time (rescaled, real: t/dstar)
c     xs      shift since t zero (0 for spatial cases)
c     nxin    collocation points in x
c     nypin   collocation points in y
c     nzcin   collocation points in z
c     nfzsin  symmetry flag (0: no symmetry)
c     fltype  flow type
c            -3: asymptotic suction boundary layer
c            -2: temporal Falkner-Skan-Cooke boundary layer
c            -1: temporal Falkner-Skan boundary layer
c             1: temporal Poiseuille flow
c             2: temporal Couette flow
c             3: temporal Blasius boundary layer
c             4: spatial Poiseuille flow
c             5: spatial Couette flow
c             6: spatial Blasius
c             7: spatial Falkner-Skan
c             8: spatial Falkner-Skan-Cooke
c             9: spatial parallel boundary layer
c     dstar   length scale (=2/yl)
c
 1010 continue

      if (fltype.eq.-3) write(*,*) 'Flowtype                      : '//
     &     'Temporal suction boundary layer (fltype = -3)'
      if (fltype.eq.-2) write(*,*) 'Flowtype                      : '//
     &     'Temporal Falkner-Skan-Cooke boundary layer (fltype = -2)'
      if (fltype.eq.-1) write(*,*)
     &     'Temporal Falkner-Skan boundary layer (fltype = -1)'
      if (fltype.eq.1) write(*,*) 'Flowtype                       : '//
     &     'Temporal Poiseuille flow (fltype =  1)'
      if (fltype.eq.2) write(*,*) 'Flowtype                       : '//
     &     'Temporal Couette flow (fltype =  2)'
      if (fltype.eq.3) write(*,*) 'Flowtype                       : '//
     &     'Temporal Blasius boundary layer (fltype =  3)'
      if (fltype.eq.4) write(*,*) 'Flowtype                       : '//
     &     'Spatial Poiseuille flow (fltype =  4)'
      if (fltype.eq.5) write(*,*) 'Flowtype                       : '//
     &     'Spatial Couette flow (fltype =  5)'
      if (fltype.eq.6) write(*,*) 'Flowtype                       : '//
     &     'Spatial Blasius boundary layer (fltype =  6)'
      if (fltype.eq.7) write(*,*) 'Flowtype                       : '//
     &     'Spatial Falkner-Skan boundary layer (fltype =  7)'
      if (fltype.eq.8) write(*,*) 'Flowtype                       : '//
     &     'Spatial Falkner-Skan-Cooke boundary layer (fltype =  8)'
      if (fltype.eq.9) write(*,*) 'Flowtype                       : '//
     &     'Spatial parallel boundary layer (fltype =  9)'

c     Check all flow types???
      if (fltype.ne.3.and.fltype.ne.-1.and.fltype.ne.-2
     &     .and.fltype.ne.1.and.fltype.ne.-3
     &     .and.(fltype.lt.6.or.fltype.gt.9)) then
         write(*,*)
         write(*,*) 'Warning! The input file contains '//
     &              'a flow type that is not verified'
         write(*,*)
c         stop
      end if
c
c     Read additional info for specific flow types
c
      rlam = 0.
      spanv = 0.
      if (fltype.eq.-1.or.fltype.eq.-2) read(12) rlam
      if (fltype.eq.6) read(12) bstart,blength
      if (fltype.ge.7) read(12) bstart,blength,rlam,spanv
      if (abs(fltype).eq.20) read(12) (gr(i),i=1,scalar)
c
c     Rescale flow data
c
      re = re * dstar
      xl = xl / dstar
      zl = zl / dstar
      yl = 2. / dstar
      t  = t  / dstar

      write(*,*) 'Re                             : ',re
      write(*,*) 't                              : ',t
      do i=1,scalar
         write(*,*) 'Pr,m1                          : ',pr(i),m1(i)
      end do
      if (abs(fltype).eq.20)  then
         do i=1,scalar
            write(*,*) 'Gr                             : ',gr(i)
         end do
      end if
      write(*,*) 'Resolution (file)              : ',nxin,nypin,nzin
      write(*,'(A,3F12.5)') ' Box                            : ',
     &     xl,yl,zl
c
c     Change internal resolution
c     
      nx=nxin*xfac
      nyp=nypin
      nz=nzin*zfac
c
c     Set bounds to default if not subset
c
      if (.not. subset) then
         bounds(1) = 1
         bounds(2) = nx
         bounds(3) = 1
         bounds(4) = nyp
         bounds(5) = 1
         bounds(6) = nz
      end if

      if (fltype.ne.1.and.ymin.lt.0.0) then
         write(*,*) 'Negative ymin given for BL, setting to zero!'
         ymin=0.0
      end if

      if (.not. uniform_y) then
         write(*,*) 'Resolution (post)              : ',nx,nyp,nz
      else
         write(*,*) 'Uniform grid in y using Clenshaw interpolation'
         write(*,*) 'Resolution (post)              : ',nx,nypc,nz
         write(*,*) 'Box y-bounds (post)            : ',ymin,'-',ymax 
         bounds(4) = nypc
      end if
c
c     Determine if vorticity is to be computed
c
      if (vort) then
         ivort = 1
      else
         ivort = 0
      end if
c
c     Allocate data structure
c
      allocate(ur(nx/2+1,nyp,nz,3+scalar),ui(nx/2+1,nyp,nz,3+scalar))

      write(*,*)
      write(*,*) 'Read data...'
c
c     Read data
c     Note that the wall is at index nyp
c
      allocate(urx(nx))
      do ll=1,3+scalar
         write(*,*) '  - component ',ll,' of ',3+scalar
         do k=1,nz
            do j=1,nyp

               urx = 0.0

               if (k.gt.(nzin+1)/2.and.k.le.nz-nzin/2) then
               else
                  read(12) (urx(i),i=1,nxin)
               end if

               do i=1,nx/2
                  ur(i,j,k,ll) = urx(2*i-1)
                  ui(i,j,k,ll) = urx(2*i)
               end do
               ur(nx/2+1,j,k,ll) = 0.
               ui(nx/2+1,j,k,ll) = 0.

               if (k.eq.nz-nzin/2+1) then
                  do i=1,nx/2
                     ur(i,j,k,ll) = 0.
                     ui(i,j,k,ll) = 0.
                  end do
               end if

            end do
         end do
      end do

      deallocate(urx)
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
c      beta(nz/2+1)=0.0
c
c     Do the shift of the data in Fourier space in the streamwise direction
c     so that the fringe region (if any) ends up in the end of the domain
c     and not in the center
c
      write(*,*) 'Shift data...'
      if (fltype.eq.-3) then
         xs = 0.
      else
         xs = xl/2.
      end if
      zs = 0.

c      xs = 0.
c      zs = zl/2

      allocate(cx(nx/2),sx(nx/2))
      do i=1,nx/2
         argx = -xs*alfa(i)
         cx(i) = cos(argx)
         sx(i) = sin(argx)
      end do
!$OMP PARALLEL PRIVATE(k,argz,i,ca,sa,ll,j,hr)
      allocate(ca(nx/2),sa(nx/2))
!$omp do
      do k=1,nz
         argz = -zs*beta(k)
         do i=1,nx/2
            ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz)
            sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz)
         end do
         do ll=1,3+scalar
            do j=1,nyp
               do i=1,nx/2
                  hr=ur(i,j,k,ll)*ca(i)-ui(i,j,k,ll)*sa(i)
                  ui(i,j,k,ll)=ui(i,j,k,ll)*ca(i)+ur(i,j,k,ll)*sa(i)
                  ur(i,j,k,ll)=hr
               end do
            end do
         end do
      end do
      deallocate(ca,sa)
!$omp end parallel

      deallocate(cx,sx)
c-------------------------------------------------------
c
c     Compute velocity derivatives
c
      write(*,*) 'Compute derivatives...'
      allocate(dur(nx/2+1,nyp,nz,3,3),dui(nx/2+1,nyp,nz,3,3))
c
c     Wall-parallel derivatives
c
!$OMP PARALLEL DO PRIVATE(ll,k,j,i)
      do k=1,nz
         do ll=1,3
            do j=1,nyp
               do i=1,nx/2
                  dur(i,j,k,ll,1)=-alfa(i)*ui(i,j,k,ll)
                  dui(i,j,k,ll,1)= alfa(i)*ur(i,j,k,ll)
                  dur(i,j,k,ll,3)=-beta(k)*ui(i,j,k,ll)
                  dui(i,j,k,ll,3)= beta(k)*ur(i,j,k,ll)
               end do
            end do
         end do
      end do
c
c     Compute wall-normal derivatives
c
      allocate(w3(nx/2*nz*nyp,itot))
      allocate(prey(nyp*2+15))
      call vcosti(nyp,prey,i)
c
c     Do the (normalized) Chebyshev transform
c
!$OMP PARALLEL DO PRIVATE(ll,k,j,i)
      do k=1,nz
         do ll=1,3
            do j=1,nyp
               do i=1,nx/2+1
                  ur(i,j,k,ll)=ur(i,j,k,ll)*(2./real(nyp-1))
                  ui(i,j,k,ll)=ui(i,j,k,ll)*(2./real(nyp-1))
               end do
            end do
         end do
      end do

!$OMP PARALLEL DO PRIVATE(k,ll)
      do k=1,nz
         do ll=1,3
            call vchbf(ur(1,1,k,ll),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
            call vchbf(ui(1,1,k,ll),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
         end do
      end do
c
c     Do the derivative
c
!$OMP PARALLEL DO PRIVATE(k,ll)
      do k=1,nz
         do ll=1,3
            call dcheb(dur(1,1,k,ll,2),ur(1,1,k,ll),nyp,nx/2+1,nx/2+1)
            call dcheb(dui(1,1,k,ll,2),ui(1,1,k,ll),nyp,nx/2+1,nx/2+1)
         end do
      end do
c
c     Chebyshev back transform and normalize
c
!$OMP PARALLEL DO PRIVATE(k,ll)
      do k=1,nz
         do ll=1,3
            call vchbb(ur(1,1,k,ll),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
            call vchbb(ui(1,1,k,ll),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
            call vchbb(dur(1,1,k,ll,2),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
            call vchbb(dui(1,1,k,ll,2),w3(1,omp_get_thread_num()+1),
     &           nyp,nx/2+1,nx/2+1,1,prey)
         end do
      end do

      deallocate(w3)

      write(*,*) 'Do normalisation...'

      dur(:,:,:,:,2) = dur(:,:,:,:,2)/yl*2.
      dui(:,:,:,:,2) = dui(:,:,:,:,2)/yl*2.


c
c     Do backward FFT to physical space
c
      write(*,*) 'Do inverse FFT...'
      allocate(prexn(nx+25))
      allocate(prezn(nz*2+25))
      allocate(wr((nx/2+1)*nyp*nz,itot))
      allocate(wi((nx/2+1)*nyp*nz,itot))

      call vrffti(nx,prexn,i)
      call vcffti(nz,prezn,i)

!$OMP PARALLEL DO PRIVATE(j,ll)
      do j=1,nyp
         do ll=1,3+scalar
            call vcfftb(ur(1,j,1,ll),ui(1,j,1,ll),
     &           wr(1,omp_get_thread_num()+1),
     &           wi(1,omp_get_thread_num()+1),
     &           nz,nx/2,(nx/2+1)*nyp,1,prezn)
            call vrfftb(ur(1,j,1,ll),ui(1,j,1,ll),
     &           wr(1,omp_get_thread_num()+1),
     &           wi(1,omp_get_thread_num()+1),
     &           nx,nz,1,(nx/2+1)*nyp,prexn)
         end do
      end do

!$OMP PARALLEL DO PRIVATE(j,ll,kk)
      do j=1,nyp
         do ll=1,3
            do kk=1,3
               call vcfftb(dur(1,j,1,ll,kk),dui(1,j,1,ll,kk),
     &              wr(1,omp_get_thread_num()+1),
     &              wi(1,omp_get_thread_num()+1),
     &              nz,nx/2,(nx/2+1)*nyp,1,prezn)
               call vrfftb(dur(1,j,1,ll,kk),dui(1,j,1,ll,kk),
     &              wr(1,omp_get_thread_num()+1),
     &              wi(1,omp_get_thread_num()+1),
     &              nx,nz,1,(nx/2+1)*nyp,prexn)
            end do
         end do
      end do

      deallocate(wr,wi)

c-------------------------------------------------------------

c
c     Compute wall-normal coordinates
c
      allocate(eta(nyp))
      if (fltype.eq.1.or.fltype.eq.2) then
c     Channel distribution
         do j=1,nyp
            eta(j) = cos(pi*real(j-1)/real(nyp-1))
         end do
      else if (fltype.eq.6.or.fltype.eq.-3.or.fltype.eq.3) then
c     Boundary-layer distribution
         do j=1,nyp
            eta(j) = yl-yl/2.*
     &           (1-cos(pi*real(j-1)/real(nyp-1)))
         end do
      else
         write(*,*) 'Check point distribution'
         stop
      end if
c
c     Compute streamwise coordinates
c
      allocate(xcoord(nx))
      do j=1,nx
         xcoord(j) = (j-1)*xl/nx
      end do
c
c     Compute spanwise coordinates
c
      allocate(zcoord(nz))
      do j=1,nz
         zcoord(j) = (j-1)*zl/nz
      end do

c--------------------------------------------------------

c
c     Compute derived quantities
c
      allocate(ar(nx/2+1,nyp,nz,2),ai(nx/2+1,nyp,nz,2))
c
c     Disturbance velocity
c
      write(*,*) 'Compute disturbance velocity...'

      if (fltype.eq.1) then

!$OMP PARALLEL DO PRIVATE(i,j,umeanr,umeani,k)
         do j=1,nyp
            umeanr = 0.
            do i=1,nx/2
               do k=1,nz
                  umeanr = umeanr + ur(i,j,k,1) + ui(i,j,k,1)
               end do
            end do
            umeanr = umeanr/real(nz)/real(nx)
            do i=1,nx/2
               do k=1,nz
                  ar(i,j,k,1) = ur(i,j,k,1) - umeanr
                  ai(i,j,k,1) = ui(i,j,k,1) - umeanr
               end do
            end do
         end do
         
      else

!$OMP PARALLEL DO PRIVATE(i,j,umeanr,umeani,k)
         do j=1,nyp
            do i=1,nx/2
               umeanr = 0.
               umeani = 0.
               do k=1,nz
                  umeanr = umeanr + ur(i,j,k,1)
                  umeani = umeani + ui(i,j,k,1)
               end do
               umeanr = umeanr/real(nz)
               umeani = umeani/real(nz)
               do k=1,nz
                  ar(i,j,k,1) = ur(i,j,k,1) - umeanr
                  ai(i,j,k,1) = ui(i,j,k,1) - umeani
               end do
            end do
         end do
      end if

c
c     Streak amplitude
c
c      do i=1,nx/2
c         uminr = 1000.
c         umaxr = -1000.
c         umini = 1000.
c         umaxi = -1000.
c         do j=1,nyp
c            do k=1,nz
c               uminr = min(uminr,ur(i,j,k,1))
c               umaxr = max(umaxr,ur(i,j,k,1))
c               umini = min(umini,ui(i,j,k,1))
c               umaxi = max(umaxi,ui(i,j,k,1))
c            end do
c         end do
c         write(20,*) (2*i-2)*xl/nx,(umaxr-uminr)/2
c         write(20,*) (2*i-1)*xl/nx,(umaxi-umini)/2
c      end do
c      stop

c
c     Disturbance dudy
c
c      do j=1,nyp
c         do i=1,nx/2
c            umeanr = 0.
c            umeani = 0.
c            do k=1,nz
c               umeanr = umeanr + dur(i,j,k,1,2)
c               umeani = umeani + dui(i,j,k,1,2)
c            end do
c            umeanr = umeanr/real(nz)
c            umeani = umeani/real(nz)
c            do k=1,nz
c               dur(i,j,k,1,2) = dur(i,j,k,1,2) - umeanr
c               dui(i,j,k,1,2) = dui(i,j,k,1,2) - umeani
c            end do
c         end do
c      end do
c
c     lambda_2
c

      allocate(omr(nx/2+1,nyp,nz,4*ivort))
      allocate(omi(nx/2+1,nyp,nz,4*ivort))

      write(*,*) 'Compute lambda2 and Q...'
c      call comp_lam(nx,nyp,nz,dur,ar(1,1,1,2),ar(1,1,1,3),ar(1,1,1,4))
c      call comp_lam(nx,nyp,nz,dui,ai(1,1,1,2),ai(1,1,1,3),ai(1,1,1,4))
      call comp_lam(nx,nyp,nz,dur,ar(1,1,1,2),ar(1,1,1,1),
     &              ar(1,1,1,1),omr,ivort)
      call comp_lam(nx,nyp,nz,dui,ai(1,1,1,2),ai(1,1,1,1),
     &              ai(1,1,1,1),omi,ivort)

      deallocate(dur,dui)

c     
c     Scaling of lambda2
c
      if (0.eq.1) then
         write(*,*) 'Scaling lambda2 in inner units'
         allocate(scalings(nx,8))
         open(unit=155,file='/home/x_phisc/bl/4000/l2scaling.data')
         do i=1,nx
            read(155,*) scalings(i,1:8)
         end do
         close(155)

c     5 is inner
c     6 is outer
c     7 is mixed 1
c     8 is mixed 2 (best)

         do k=1,nz
            do j=1,nyp
               do i=1,nx/2
                  ar(i,j,k,2)=ar(i,j,k,2) / scalings(i*2-1,5)
                  ai(i,j,k,2)=ai(i,j,k,2) / scalings(i*2  ,5)

c                 ur(i,j,k,1) = scalings(i*2-1,5)
c                 ui(i,j,k,1) = scalings(i*2,5)

               end do
            end do
         end do
         
         deallocate(scalings)

      end if




c
c     Compute the values on even gridpoints for volume rendering 
c     using the Clenshaw formula
c
      if (uniform_y) then

         write(*,*) "Interpolating data on a uniform grid in y..."

         allocate(urc(nx/2+1,nypc,nz,3+scalar))
         allocate(uic(nx/2+1,nypc,nz,3+scalar))
         allocate(arc(nx/2+1,nypc,nz,2),aic(nx/2+1,nypc,nz,2))
c
c     velocity
c
         if (wu.eq.1) then
            write(*,*) '   Interpolating u...'
            call clenshaw(nx,nyp,nypc,nz,ur(1,1,1,1),ui(1,1,1,1),
     &           urc(1,1,1,1),uic(1,1,1,1),yl,ymin,ymax,dy,fltype,itot)
         end if
         if (wv.eq.1) then
            write(*,*) '   Interpolating v...'
            call clenshaw(nx,nyp,nypc,nz,ur(1,1,1,2),ui(1,1,1,2),
     &           urc(1,1,1,2),uic(1,1,1,2),yl,ymin,ymax,dy,fltype,itot)
         end if
         if (ww.eq.1) then
            write(*,*) '   Interpolating w...'
            call clenshaw(nx,nyp,nypc,nz,ur(1,1,1,3),ui(1,1,1,3),
     &           urc(1,1,1,3),uic(1,1,1,3),yl,ymin,ymax,dy,fltype,itot)
         end if
         if (ws.eq.1) then
            write(*,*) '   Interpolating t...'
            call clenshaw(nx,nyp,nypc,nz,ur(1,1,1,4),ui(1,1,1,4),
     &           urc(1,1,1,4),uic(1,1,1,4),yl,ymin,ymax,dy,fltype,itot)
         end if
c
c     lambda2
c
         if (wl2.eq.1) then
            write(*,*) '   Interpolating lambda2...'
            call clenshaw(nx,nyp,nypc,nz,ar(1,1,1,2),ai(1,1,1,2),
     &           arc(1,1,1,2),aic(1,1,1,2),yl,ymin,ymax,dy,fltype,itot)         
         end if
c
c     ufluct
c
         if (wf.eq.1) then
            write(*,*) '   Interpolating ufluct...'
            call clenshaw(nx,nyp,nypc,nz,ar(1,1,1,1),ai(1,1,1,1),
     &           arc(1,1,1,1),aic(1,1,1,1),yl,ymin,ymax,dy,fltype,itot)         
         end if
c
c     vorticity
c 
         if (vort) then

            allocate(omrc(nx/2+1,nypc,nz,3))
            allocate(omic(nx/2+1,nypc,nz,3))

            if (wo1.eq.1) then
               write(*,*) '   Interpolating omega_x...'
               call clenshaw(nx,nyp,nypc,nz,omr(1,1,1,1),omi(1,1,1,1),
     &              omrc(1,1,1,1),omic(1,1,1,1),yl,ymin,ymax,dy,
     &              fltype,itot)
            end if
            if (wo2.eq.1) then
               write(*,*) '   Interpolating omega_y...'
               call clenshaw(nx,nyp,nypc,nz,omr(1,1,1,2),omi(1,1,1,2),
     &              omrc(1,1,1,2),omic(1,1,1,2),yl,ymin,ymax,dy,
     &              fltype,itot)
            end if
            if (wo3.eq.1) then
               write(*,*) '   Interpolating omega_z...'
               call clenshaw(nx,nyp,nypc,nz,omr(1,1,1,3),omi(1,1,1,3),
     &              omrc(1,1,1,3),omic(1,1,1,3),yl,ymin,ymax,dy,
     &              fltype,itot)
            end if
         end if
c
c     scalars
c
         do i=1,numscalars
            write(*,*) '   Interpolating scalar',scalars(i),'...'
            call clenshaw(nx,nyp,nypc,nz,ur(1,1,1,3+scalars(i)),
     &           ui(1,1,1,3+scalars(i)),urc(1,1,1,3+scalars(i)),
     &           uic(1,1,1,3+scalars(i)),yl,ymin,ymax,dy,fltype,itot)
         end do
      end if      
c
c     Write data for different postprocessing softwares
c

      write(*,*)
      write(*,*) 'Output format:',type

      if (type.eq.'dx') then
         write(*,*) 'dx: OpenDX'
         call write_dx(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xl,zl,
     &        t,fileout,scalar)
      elseif (type.eq.'vt') then
         write(*,*) 'vt: vtk (paraview, visit)'
         write(*,*) 'file name:',fileout
         if (.not. uniform_y) then 
            nypc=nyp 
            call write_xml(ur,ui,ar,ai,dur,dui,nx,nypc,nz,eta,xl,zl,
     &        t,filein,fileout,scalar,omr,omi,ivort,dy,uniform_y,
     &        bounds,vars,numvar,ymin,numscalars,scalars)
         else
            call write_vtk(urc,uic,arc,aic,dur,dui,nx,nypc,nz,eta,xl,zl,
     &        t,filein,fileout,scalar,omrc,omic,ivort,dy,uniform_y,
     &        bounds,vars,numvar,ymin,numscalars,scalars)            
         end if

      elseif (type.eq.'ef') then
         write(*,*) 'ef: EnSight Gold (Fortran binary) eg EnSight'
         cb=.false.
         call write_engold(ur,ui,dur,dui,ar,ai,nx,nyp,nz,eta,xl,yl,zl,
     &        re,rlam,spanv,fltype,t,xcoord,zcoord,fileout,cb)
      elseif (type.eq.'ec') then
         write(*,*) 'ec: EnSight Gold (C binary), eg Paraview, EnSight'
         cb=.true.
         call write_engold(ur,ui,dur,dui,ar,ai,nx,nyp,nz,eta,xl,yl,zl,
     &        re,rlam,spanv,fltype,t,xcoord,zcoord,fileout,cb)
#ifdef HDF
      elseif (type.eq.'df') then
         write(*,*) 'df: HDF'
         call write_hdf(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xcoord,zcoord,
     &        t,fileout,scalar)
#endif
      else
         write(*,*) 'Not a valid output option.'
         stop
      end if

c      deallocate(dur,dui)
      deallocate(ar,ai)
      deallocate(ur,ui)
      if (vort) deallocate(omr,omi)

      if (uniform_y) then
         deallocate(urc,uic)
         deallocate(arc,aic)
         if (vort) then
            deallocate(omrc,omic)
         end if
      end if

      end program lambda2
