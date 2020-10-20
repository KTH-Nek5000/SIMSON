c ***********************************************************************
C
c $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/lambda2/write_vtk.f $
c $LastChangedDate: 2011-05-26 14:29:53 +0200 (Thu, 26 May 2011) $
c $LastChangedBy: ilak@MECH.KTH.SE $
c $LastChangedRevision: 1671 $
c
c ***********************************************************************
      subroutine write_xml(ur,ui,ar,ai,dur,dui,nx,nyp,nz,eta,xl,zl,
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
      integer wu,wv,ww,wl2,wo1,wo2,wo3,wf
      integer scalars(numscalars)
      integer,parameter :: sfloat=4
      
      integer bytes_x,bytes_y,bytes_z,bytes_f
      integer offset

      if (uniform_y) then
         write(*,*) 'Uniform grid for XML not yet done.'
         stop
      end if

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
c      open(unit=22,file=trim(fileout)//'.vtr',form='binary')
c      open(unit=22,file=trim(fileout)//'.vtr',form='binary',
c     &     convert='little_endian')
      open(unit=22,file=trim(fileout)//'.vtr',form='unformatted',
     &     convert='little_endian',access='stream')

      nn = char(10)

      write(22) '<!-- This is a SIMSON file -->'//nn
      write(22) '<VTKFile type="RectilinearGrid" version="0.1" '//
     &     'byte_order="LittleEndian">'//nn


      write(string,'(a,6i5,a)') 
     &     '  <RectilinearGrid WholeExtent="',
     &     0,(ie-is+1)*2-1,0,je-js,0,ze-zs,'">'//nn
      write(22) trim(string)

      write(22) '    <FieldData>'//nn
      write(22) '      <DataArray type="Float32" Name="TIME" '//
     &     'NumberOfTuples="1" format="ascii">'//nn
      write(string,'(e18.7,a)') t,nn
      write(22) trim(string)
      write(22) '      </DataArray>'//nn
      write(22) '    </FieldData>'//nn

      write(string,'(a,6i5,a)') 
     &     '    <Piece Extent="',
     &     0,(ie-is+1)*2-1,0,je-js,0,ze-zs,'">'//nn
      write(22) trim(string)

      write(22) '      <Coordinates>'//nn

      offset = 0
      bytes_f = ((ze-zs)+1) * ((je-js)+1) * (ie-is+1)*2 * sfloat
c
c     x-coordinates
c
      write(string,'(a,i10,a)') '        <DataArray type="Float32" '//
     &     'Name="X_COORDINATES" format="appended" offset="',
     &     offset,'"/>'//nn
      write(22) trim(string)
      bytes_x = (ie-is+1)*2 * sfloat
      offset = offset + bytes_x + 4
c
c     y-coordinates
c
      write(string,'(a,i10,a)') '        <DataArray type="Float32" '//
     &     'Name="Y_COORDINATES" format="appended" offset="',
     &     offset,'"/>'//nn
      write(22) trim(string)
      bytes_y = ((je-js)+1) * sfloat
      offset = offset + bytes_y + 4
c
c     z-coordinates
c
      write(string,'(a,i10,a)') '        <DataArray type="Float32" '//
     &     'Name="Z_COORDINATES" format="appended" offset="',
     &     offset,'"/>'//nn
      write(22) trim(string)
      bytes_z = ((ze-zs)+1) * sfloat
      offset = offset + bytes_z + 4

      write(22) '      </Coordinates>'//nn

      write(22) '      <PointData Scalars="vel_u">'//nn

      if (wu.eq.1) then
c
c     u-velocity
c
         write(*,*) 'Writing u velocity...'
         write(string,'(a,i10,a)') '        <DataArray '//
     &        'type="Float32" '//
     &        'Name="vel_u" format="appended" offset="',
     &        offset,'"/>'//nn
         write(22) trim(string)
         offset = offset + bytes_f + 4
      end if
      if (wv.eq.1) then
c
c     v-velocity
c
         write(*,*) 'Writing v velocity...'
         write(string,'(a,i10,a)') '        <DataArray '//
     &        'type="Float32" '//
     &        'Name="vel_v" format="appended" offset="',
     &        offset,'"/>'//nn
         write(22) trim(string)
         offset = offset + bytes_f + 4
      end if
      if (ww.eq.1) then
c
c     w-velocity
c
         write(*,*) 'Writing w velocity...'
         write(string,'(a,i10,a)') '        <DataArray '//
     &        'type="Float32" '//
     &        'Name="vel_w" format="appended" offset="',
     &        offset,'"/>'//nn
         write(22) trim(string)
         offset = offset + bytes_f + 4
      end if
      if (wf.eq.1) then
c
c     ufluct
c
         write(*,*) 'Writing ufluct velocity...'
         write(string,'(a,i10,a)') '        <DataArray '//
     &        'type="Float32" '//
     &        'Name="ufluct" format="appended" offset="',
     &        offset,'"/>'//nn
         write(22) trim(string)
         offset = offset + bytes_f + 4
      end if
      if (wl2.eq.1) then
c
c     lambda2
c
         write(*,*) 'Writing lambda2 velocity...'
         write(string,'(a,i10,a)') '        <DataArray '//
     &        'type="Float32" '//
     &        'Name="lambda2" format="appended" offset="',
     &        offset,'"/>'//nn
         write(22) trim(string)
         offset = offset + bytes_f + 4
      end if
      if (ivort.gt.0) then
c
c     vorticity 1
c
         if (wo1.eq.1) then
            write(*,*) 'Writing vort_x velocity...'
            write(string,'(a,i10,a)') '        <DataArray '//
     &           'type="Float32" '//
     &           'Name="vort_x" format="appended" offset="',
     &           offset,'"/>'//nn
            write(22) trim(string)
            offset = offset + bytes_f + 4
         end if
c
c     vorticity 2
c
         if (wo2.eq.1) then
            write(*,*) 'Writing vort_y velocity...'
            write(string,'(a,i10,a)') '        <DataArray '//
     &           'type="Float32" '//
     &           'Name="vort_y" format="appended" offset="',
     &           offset,'"/>'//nn
            write(22) trim(string)
            offset = offset + bytes_f + 4
         end if
c
c     vorticity 3
c
         if (wo3.eq.1) then
            write(*,*) 'Writing vort_z velocity...'
            write(string,'(a,i10,a)') '        <DataArray '//
     &           'type="Float32" '//
     &           'Name="vort_z" format="appended" offset="',
     &           offset,'"/>'//nn
            write(22) trim(string)
            offset = offset + bytes_f + 4
         end if
      end if


      if (scalar.gt.0) then
c
c     scalars
c
         do ii=1,numscalars
            write(*,*) 'Writing scalar number ',scalars(ii),'...'
            write(string,'(a,i1,a,i10,a)') '        <DataArray '//
     &           'type="Float32" '//
     &           'Name="scal',
     &           scalars(ii),
     &           '" format="appended" offset="',
     &           offset,'"/>'//nn
            write(22) trim(string)
            offset = offset + bytes_f + 4
         end do
      end if

c
c     Footer
c
      write(22) '      </PointData>'//nn
      write(22) '      <CellData>'//nn
      write(22) '      </CellData>'//nn
      write(22) '    </Piece>'//nn
      write(22) '  </RectilinearGrid>'//nn
c
c     Data
c
      write(22) '  <AppendedData encoding="raw">'//nn
      write(22) '_'
      write(22) bytes_x
      write(22) (real(xl/nx*(2*i-2),4),real(xl/nx*(2*i-1),4),
     &     i=is,ie)
      write(22) bytes_y
      write(22) (real(eta(j),4),j=js,je)
      write(22) bytes_z
      write(22) (real(-zl/2+zl/nz*(k-1),4),k=zs,ze)

      if (wu.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(ur(i,j,k,1),4),real(ui(i,j,k,1),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wv.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(ur(i,j,k,2),4),real(ui(i,j,k,2),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (ww.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(ur(i,j,k,3),4),real(ui(i,j,k,3),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wf.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(ar(i,j,k,1),4),real(ai(i,j,k,1),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wl2.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(ar(i,j,k,2),4),real(ai(i,j,k,2),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wo1.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(omr(i,j,k,1),4),real(omi(i,j,k,1),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wo2.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(omr(i,j,k,2),4),real(omi(i,j,k,2),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (wo3.eq.1) then
         write(22) bytes_f
         write(22)
     &        ((( real(omr(i,j,k,3),4),real(omi(i,j,k,3),4),
     &        i=is,ie),j=js,je),k=zs,ze)
      end if
      if (scalar.gt.0) then
         do ii=1,numscalars
            write(22) bytes_f
            write(22)
     &           ((( real(ur(i,j,k,3+scalars(ii)),4),
     &           real(ui(i,j,k,3+scalars(ii)),4),
     &           i=is,ie),j=js,je),k=zs,ze)
         end do
      end if


      write(22) nn


      write(22) '  </AppendedData>'//nn
      write(22) '</VTKFile>'//nn
      close(22)

      end subroutine write_xml
