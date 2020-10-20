c ***********************************************************************
C
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine write_vtk(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xl,zl,
     &     t,filein,fileout,scalar,omr,omi,ivort,dy,uniform_y,
     &     bounds,vars,numvar,ymin,numscalars,scalars)
c
c     Creates OpenDX readable files
c
      implicit none

      integer i,j,k,nx,nyp,nz,ii
      integer is,ie,js,je,scalar,ivort,zs,ze
      real*8 xl,zl,t
      real ur(nx/2+1,nyp,nz,3+scalar),ui(nx/2+1,nyp,nz,3+scalar)
      real ar(nx/2+1,nyp,nz,2),ai(nx/2+1,nyp,nz,2)
      real omr(nx/2+1,nyp,nz,3*ivort),omi(nx/2+1,nyp,nz,3*ivort)
      real dur(nx/2+1,nyp,nz,3,3),dui(nx/2+1,nyp,nz,3,3)
      real eta(nyp)
      real dx,dy,dz,zhalf,ymin
      character*80 filein,fileout
      character*1000 string
      character(len=100) path_name
      character(len=20) t_string
      character nn
      logical uniform_y,fltype
      integer bounds(6)
      integer numvar,numscalars
      integer vars(numvar)
      integer wu,wv,ww,wl2,wo1,wo2,wo3,wo4,wf
      integer scalars(numscalars),scwrite(numscalars)
c
c     figure out which variables to write
c
      wu=0
      wv=0
      ww=0
      wf=0
      wl2=0
      wo1=0
      wo2=0
      wo3=0
      do i=1,numvar
         if (vars(i) .eq. 1) then
            wu=1
         elseif (vars(i) .eq. 2) then
            wv=1
         elseif (vars(i) .eq. 3) then
            ww=1
         elseif (vars(i) .eq. 4) then
            wl2=1
         elseif (vars(i) .eq. 5) then
            wo1=1
         elseif (vars(i) .eq. 6) then
            wo2=1
         elseif (vars(i) .eq. 7) then
            wo3=1
         elseif (vars(i) .eq. 8) then
            wf=1
         elseif (vars(i) .eq. 9) then
            wo4=1
         end if
      end do
c
c     Complete domain
c
      is = bounds(1)
      ie = bounds(2)

      js = bounds(3)
      je = bounds(4)

      zs = bounds(5)
      ze = bounds(6)
c
c     Some restrictions in output if desired
c     Note that j=1 is the top boundary, j=nyp the lower wall
c     is should be odd, ie even
c      is = 3501
c      ie = 4500

      write(*,*) 'x: ',is,' - ',ie,'(1 -',nx,')'
      write(*,*) 'y: ',js,' - ',je,'(1 -',nyp,')'
      write(*,*) 'z: ',zs,' - ',ze,'(1 -',nz,')'

      is = (is-1)/2+1
      ie = ie/2
c
c     VTK Format
c
c      open(unit=22,file=trim(fileout)//'.vtk',form='binary',
c     &     convert='big_endian')
      open(unit=22,file=trim(fileout)//'.vtk',form='unformatted',
     &     convert='big_endian',ACCESS='stream')

      nn = char(10)

      write(22) "# vtk DataFile Version 2.0"//nn
      write(string,'(a,f25.10,a)') 
     &     "SIMSON velocity field t=",t,nn
      write(22) trim(string)
      write(22) "BINARY"//nn

      if (.not. uniform_y) then
         write(22) "DATASET RECTILINEAR_GRID"//nn
      else
         write(22) "DATASET STRUCTURED_POINTS"//nn
      end if

      write(string,'(a,3i5,a)') 
     &     "DIMENSIONS ",(ie-is+1)*2,je-js+1,(ze-zs+1),nn
      write(22) trim(string)

      if (.not. uniform_y) then

         write(string,'(a,i5,a,a)') 
     &        "X_COORDINATES ",(ie-is+1)*2," float",nn
         write(22) trim(string)
         write(22) (real(xl/nx*(2*i-2),4),
     &        real(xl/nx*(2*i-1),4),
     &        i=is,ie)
         
         write(string,'(a,i5,a,a)') 
     &        "Y_COORDINATES ",je-js+1," float",nn
         write(22) trim(string)
         write(22) (real(eta(j),4),j=js,je)

         write(string,'(a,i5,a,a)')
     &        "Z_COORDINATES ",(ze-zs+1)," float",nn
         write(22) trim(string)
         write(22) (real(-zl/2+zl/nz*(k-1),4),
     &        k=zs,ze)
         
      else

         zhalf = -zl/2
         write(string,'(a,f30.15,f30.15,f30.15,a)')
     &     "ORIGIN ",real(xl/nx*(2*is-2),4), 
     &               real(ymin+dy*(js-1),4),
     &               real(zhalf+zl/nz*(zs-1),4),nn   
         write(22) trim(string)

         dx = xl/real(nx)
         dz = zl/real(nz)

         write(*,*) "Spacing of uniform grid:"
         write(*,*) "dx:",dx ,"dy:",dy ,"dz:",dz      
         write(string,'(a,3f30.15,a)')
     &     "SPACING ",real(dx,4),real(dy,4),real(dz,4),nn
c      write(*,*),string

         write(22) trim(string)

      end if

      write(22) "FIELD MetaData 4"//nn

      write(22) "TIME 1 1 float"//nn
      write(22) real(t,4)

      write(22) "FILE 1 100 char"//nn
      write(path_name,'(a)') trim(filein)
      write(22) path_name

      write(22) "PATH 1 100 char"//nn
      call getenv('PWD',path_name)
      write(22) path_name

      write(22) "GENERATED 1 20 char"//nn
      call time_string(t_string)
c      t_string = '11111111111111111111'
      write(22) t_string


      write(string,'(a,i10,a)') 
     &     "POINT_DATA ",((ie-is+1)*2)*(je-js+1)*(ze-zs+1),nn
      write(22) trim(string)

      if (wu.eq.1) then
         write(*,*) "Writing u velocity..."
         write(22) "SCALARS vel_u float 1"//nn
         write(22) "LOOKUP_TABLE default"//nn
         write(22)
     &        ((( real(ur(i,j,k,1),4),real(ui(i,j,k,1),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if

      if (wv.eq.1) then
         write(*,*) "Writing v velocity..."
         write(22) "SCALARS vel_v float 1"//nn
         write(22) "LOOKUP_TABLE default"//nn
         write(22)
     &        ((( real(ur(i,j,k,2),4),real(ui(i,j,k,2),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if

      if (ww.eq.1) then
         write(*,*) "Writing w velocity..."
         write(22) "SCALARS vel_w float 1"//nn
         write(22) "LOOKUP_TABLE default"//nn
         write(22)
     &        ((( real(ur(i,j,k,3),4),real(ui(i,j,k,3),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if

      if (wf.eq.1) then
         write(*,*) "Writing ufluct..."
         write(22) "SCALARS ufluct float 1"//nn
         write(22) "LOOKUP_TABLE default"//nn
         write(22)
     &        ((( real(ar(i,j,k,1),4),real(ai(i,j,k,1),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if

      if (wl2.eq.1) then
         write(*,*) "Writing lambda2..."
         write(22) "SCALARS lambda2 float 1"//nn
         write(22) "LOOKUP_TABLE default"//nn
         write(22)
     &     ((( real(ar(i,j,k,2),4),real(ai(i,j,k,2),4),
     &     i=is,ie),j=js,je),k=zs,ze)
      end if

c
c     write out the vorticity (optional)
c
      if (ivort.gt.0) then

         if (wo1.eq.1) then
            write(*,*) "Writing vort_x..."
            write(22) "SCALARS vort_x float 1"//nn
            write(22) "LOOKUP_TABLE default"//nn
            write(22)
     &           ((( real(omr(i,j,k,1),4),real(omi(i,j,k,1),4),
     &           i=is,ie),j=js,je),k=zs,ze)
         end if   

         if (wo2.eq.1) then
            write(*,*) "Writing vort_y..."
            write(22) "SCALARS vort_y float 1"//nn
            write(22) "LOOKUP_TABLE default"//nn
            write(22)
     &           ((( real(omr(i,j,k,2),4),real(omi(i,j,k,2),4),
     &           i=is,ie),j=js,je),k=zs,ze)
         end if

         if (wo3.eq.1) then
            write(*,*) "Writing vort_z..."
            write(22) "SCALARS vort_z float 1"//nn
            write(22) "LOOKUP_TABLE default"//nn
            write(22)
     &           ((( real(omr(i,j,k,3),4),real(omi(i,j,k,3),4),
     &           i=is,ie),j=js,je),k=zs,ze)
         end if

         if (wo4.eq.1) then
            write(*,*) "Writing omega..."
            write(22) "SCALARS vort float 1"//nn
            write(22) "LOOKUP_TABLE default"//nn
            write(22)
     &           ((( real(omr(i,j,k,4),4),real(omi(i,j,k,4),4),
     &           i=is,ie),j=js,je),k=zs,ze)
         end if
      end if
c
c     write out scalar
c
      if (scalar.gt.0) then
         do ii=1,numscalars
            write(*,*) "Writing scalar number",scalars(ii),"..."
            
            write(string,'(a,i1,a,a)') 
     &           "SCALARS scalar",scalars(ii)," float 1",nn
            write(22) trim(string)
            write(22) "LOOKUP_TABLE default"//nn
            write(22)
     &           ((( real(ur(i,j,k,3+scalars(ii)),4),
     &               real(ui(i,j,k,3+scalars(ii)),4),
     &           i=is,ie),j=js,je),k=zs,ze)         
         end do
      end if
c      write(22) "SCALARS Q float 1"//nn
c      write(22) "LOOKUP_TABLE default"//nn
c      write(22)
c     &     ((( real(ar(i,j,k,3),4),real(ai(i,j,k,3),4),
c     &     i=is,ie),j=js,je),k=zs,ze)

c      write(22) "SCALARS enstrophy float 1"//nn
c      write(22) "LOOKUP_TABLE default"//nn
c      write(22)
c     &     ((( real(ar(i,j,k,4),4),real(ai(i,j,k,4),4),
c     &     i=is,ie),j=js,je),k=zs,ze)

c      write(22) "SCALARS dudy float 1"//nn
c      write(22) "LOOKUP_TABLE default"//nn
c      write(22)
c     &     ((( real(dur(i,j,k,1,2),4),real(dui(i,j,k,1,2),4),
c     &     i=is,ie),j=js,je),k=zs,ze)

      close(22)

      end subroutine write_vtk
