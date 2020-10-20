c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine write_engold(ur,ui,dur,dui,ar,ai,nx,nyp,nz,eta,xl,yl,
     &     zl,re,rlam,spanv,fltype,t,xcoord,zcoord,fileout,cb)
c
c     Creates an EnSight Gold case file
c
      integer timeset,nsteps,nstart,ninc
      integer i,j,k,ll,mm,fltype,nx,nyp,nz,reclength
      real xl,yl,zl,re,rlam,spanv,t
      real ur(nx/2+1,nyp,nz,3),ui(nx/2+1,nyp,nz,3)
      real ar(nx/2+1,nyp,nz,2),ai(nx/2+1,nyp,nz,2)
      real dur(nx/2+1,nyp,nz,3,3),dui(nx/2+1,nyp,nz,3,3)
      real xcoord(nx),eta(nyp),zcoord(nz)
      logical grad,l2,vort,cb
      character*80 fileout
      character*80 tbinary,tnodeid,telementid,textent
      character*80 trectilinear,tpart,tblock
      character*80 tdescription1,tdescription2,tdescription3

      l2=.true.
c     Vorticity not implemented yet???
      vort=.false.
      grad=.false.
c
c     Only one time step at the moment
c
      timeset=1
      nsteps=1
      nstart=1
      ninc=1
      open(unit=22,file=trim(fileout)//'.case',form='formatted')
c
c     Write format section
c
      write(22,'(A)') 'FORMAT'
      write(22,'(A)') 'type: ensight gold'
      write(22,'(A)')
c
c     Write geometry section
c
      write(22,'(A)') 'GEOMETRY'
      write(22,'(A)') 'model:    1  '//trim(fileout)//'.geo'
      write(22,'(A)')
c
c     Write variable section
c
      write(22,'(A)') 'VARIABLE'
      write(22,'(A,F8.2)') 'constant per case:  1 X-length   ',xl
      write(22,'(A,F8.2)') 'constant per case:  1 Y-length   ',yl
      write(22,'(A,F8.2)') 'constant per case:  1 Z-length   ',zl
      write(22,'(A,F8.2)') 'constant per case:  1 Re         ',re
      write(22,'(A,I8)') 'constant per case:  1 Nx         ',nx
      write(22,'(A,I8)') 'constant per case:  1 Ny         ',nyp
      write(22,'(A,I8)') 'constant per case:  1 Nz         ',nz
      write(22,'(A,I3)') 'constant per case:  1 Flowtype   ',fltype
      write(22,'(A,F8.2)') 'constant per case:  1 Rlam       ',rlam
      write(22,'(A,F8.2)') 'constant per case:  1 Spanv      ',spanv
      write(22,'(A)') 'vector per node:    1 Velocity   '//
     &     trim(fileout)//'.vel'
      if ( grad ) then
         write(22,'(A)') 'tensor asym per node:    1 Gradient   '//
     &        trim(fileout)//'.grad'
      end if
      if ( l2 ) then
         write(22,'(A)') 'scalar per node:    1 Lambda2   '//
     &        trim(fileout)//'.l2'
      end if
      if ( vort ) then
         write(22,'(A)') 'vector per node:    1 Vorticity   '//
     &        trim(fileout)//'.vort'
      end if
c
c     Write time section
c
      write(22,'(A)')
      write(22,'(A)') 'TIME'
      write(22,'(A,I5)') 'time set:                     ',timeset
      write(22,'(A,I5)') 'number of steps:              ',nsteps
      write(22,'(A,I5)') 'filename start number:        ',nstart
      write(22,'(A,I5)') 'filename increment:           ',ninc
      write(22,'(A,F8.2)') 'time values:                  ',t
c      do j=2,nsteps
c         write(22,'(F8.2)') t(j))
c         if ( mod(j,6) .eq. 0 ) then
c            write(22,'(A)') '')
c         end if
c      end do
c
c     Close the file
c
      close(22)
c
c     Write geometry and data files
c
      if (cb) then
         write(*,*) 'Generating EnSight Gold C binary files'
         tbinary='C Binary'
      else
         write(*,*) 'Generating EnSight Gold Fortran binary files'
         tbinary='Fortran Binary'
      end if
      tnodeid='node id assign'
      telementid='element id assign'
      textent='extents'
      tpart='part'
      tblock='block'
      trectilinear='block rectilinear'

      tdescription1='Structured DNS grid'
c      tdescription2=tflowtype
      tdescription2='...'
      tdescription3='Main region'

      if (cb) then
c
c     Note that one record seems to be counted as four bytes, thus divide by 4
c
         reclength=(9*80+4*4+(6+nx+nyp+nz)*4)/4
         open(unit=22,file=trim(fileout)//'.geo',form='unformatted',
     &        access='direct',recl=reclength)
         write(22,rec=1) tbinary,
     &        tdescription1,
     &        tdescription2,
     &        tnodeid,
     &        telementid,
     &        textent,
     &        real(0.,4),real(xl,4),
     &        real(0.,4),real(yl,4),
     &        real(0.,4),real(zl,4),
     &        tpart,
     &        1,
     &        tdescription3,
     &        trectilinear,
     &        nx,nyp,nz,
     &        (real(xcoord(i),4),i=1,nx),
     &        (real(eta(i),4),i=1,nyp),
     &        (real(zcoord(i),4),i=1,nz)
      else
         open(unit=22,file=trim(fileout)//'.geo',form='unformatted')

         write(22) tbinary
         write(22) tdescription1
         write(22) tdescription2
         write(22) tnodeid
         write(22) telementid
         write(22) textent
         write(22) real(0.,4),real(xl,4),
     &             real(0.,4),real(yl,4),
     &             real(0.,4),real(zl,4)
         write(22) tpart
         write(22) 1
         write(22) tdescription3
         write(22) trectilinear
         write(22) nx,nyp,nz
         write(22) (real(xcoord(i),4),i=1,nx)
         write(22) (real(eta(i),4),i=1,nyp)
         write(22) (real(zcoord(i),4),i=1,nz)
      end if
      close(22)
c
c     Write velocity file
c
      tdescription1='Velocity'
      if (cb) then
         reclength=(3*80+1*4+(3*nx*nyp*nz)*4)/4
         open(unit=22,file=trim(fileout)//'.vel',form='unformatted',
     &        access='direct',recl=reclength)
         write(22,rec=1) tdescription1,
     &        tpart,
     &        1,
     &        tblock,
     &        ((((real(ur(i,j,k,ll),4),real(ui(i,j,k,ll),4),
     &        i=1,nx/2),j=1,nyp),k=1,nz),ll=1,3)
      else
         open(unit=22,file=trim(fileout)//'.vel',form='unformatted')
         write(22) tdescription1
         write(22) tpart
         write(22) 1
         write(22) tblock
         do ll=1,3
            write(22) (((real(ur(i,j,k,ll),4),real(ui(i,j,k,ll),4),
     &           i=1,nx/2),j=1,nyp),k=1,nz)
         end do
      end if
      close(22)
c
c     Write velocity gradient tensor file
c
      if (grad) then
         tdescription1='Gradient'
         if (cb) then
            reclength=(3*80+1*4+4*(9*nx*nyp*nz))/4
            open(unit=22,file=trim(fileout)//'.grad',form='unformatted',
     &           access='direct',recl=reclength)
            write(22,rec=1) tdescription1,
     &           tpart,
     &           1,
     &           tblock,
     &           (((((real(dur(i,j,k,ll,mm),4),real(dui(i,j,k,ll,mm),4),
     &           i=1,nx/2),j=1,nyp),k=1,nz),ll=1,3),mm=1,3)
         else
            open(unit=22,file=trim(fileout)//'.grad',form='unformatted')
            write(22) tdescription1
            write(22) tpart
            write(22) 1
            write(22) tblock
            do mm=1,3
               do ll=1,3
                  write(22) (((real(dur(i,j,k,ll,mm),4),
     &                         real(dui(i,j,k,ll,mm),4),
     &                         i=1,nx/2),j=1,nyp),k=1,nz)
               end do
            end do
         end if
         close(22)
      end if
c
c     Write vorticity file
c
      if (vort) then
         tdescription1='Vorticity'
         if (cb) then
            reclength=(3*80+1*4+4*(3*nx*nyp*nz))/4
            open(unit=22,file=trim(fileout)//'.vort',form='unformatted',
     &           access='direct',recl=reclength)
            write(22,rec=1) tdescription1,
     &           tpart,
     &           1,
     &           tblock,
     &           ((((real(ur(i,j,k,ll),4),real(ui(i,j,k,ll),4),
     &           i=1,nx/2),j=1,nyp),k=1,nz),ll=1,3)
         else
            open(unit=22,file=trim(fileout)//'.vort',form='unformatted')
            write(22) tdescription1
            write(22) tpart
            write(22) 1
            write(22) tblock
            do ll=1,3
               write(22) (((real(ur(i,j,k,ll),4),real(ui(i,j,k,ll),4),
     &              i=1,nx/2),j=1,nyp),k=1,nz)
            end do
         end if
         close(22)
      end if
c
c     Write lambda2 file
c
      tdescription1='Lambda2'
      if ( l2 ) then
         if (cb) then
            reclength=(3*80+1*4+(1*nx*nyp*nz)*4)/4
            open(unit=22,file=trim(fileout)//'.l2',form='unformatted',
     &           access='direct',recl=reclength)
            write(22,rec=1) tdescription1,
     &           tpart,
     &           1,
     &           tblock,
     &           (((real(ar(i,j,k,2),4),real(ai(i,j,k,2),4),
     &           i=1,nx/2),j=1,nyp),k=1,nz)
         else
            open(unit=22,file=trim(fileout)//'.l2',form='unformatted')
            write(22) tdescription1
            write(22) tpart
            write(22) 1
            write(22) tblock
            write(22) (((real(ar(i,j,k,2),4),real(ai(i,j,k,2),4),
     &           i=1,nx/2),j=1,nyp),k=1,nz)
         end if
         close(22)
      end if

      end subroutine write_engold
