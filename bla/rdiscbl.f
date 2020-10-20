c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdiscbl(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,
     &     pr,xl,zl,t,xs,dstar,fltype,bstart,blength,rlam,m1,spanv,
     &     spat,namnin,varsiz,m,boxr,boxi,w3,urx,nxtmp,my_node,gr)
c
c     Reads m variables from file namnin and puts then into ur
c
      implicit none

      include 'par.f'

      character*80 namnin
      integer m,fltype,nxtmp
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real u0low,u0upp,w0low,w0upp,du0upp,m1(scalar),gr(scalar)
      real re,pr(scalar),xl,zl,t,xs,dstar,bstart,blength,spanv,rlam
      logical pou,varsiz,spat
      real urx(nx),urz(nz)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real w3(nx/2,mbz,nyp)

      integer x,y,z,i,zb,zs,nxz,nxr,zu
      integer nxin,nypin,nzcin,nfzsin
      real preyin(2*nyp+15),prey(2*nyp+15)

      integer my_node,zb_t,node_t,zbb
c
c
c
      if (index(namnin,'NONE').eq.1.and.len_trim(namnin).eq.4) then
               re = 16403.
               pou = .false.
               xl = 6.2831853070000
               zl = 3.141592654000
               t = 0.
               xs = 0.
               
               nxin = nx
               nypin = nyp
               nzcin = nz
               nfzsin = 0
               
               fltype= 1
               rlam  = 0.
               dstar = 1.
               spanv = 0.
               
               ur = 0.
               ui = 0.
               return
               
      end if
c     Open input file, read file header
c
      open(unit=12,file=namnin,status='old',form='unformatted')
c
c     Description of parameters:
c     re      Reynolds number (rescaled, real: re*dstar)
c     pr      Prandtl number (only read for scalar=1)
c     pou     Not used anymore, always set to false
c     xl      streamwise length (rescaled, real: xl/dstar)
c     zl      spanwise width (rescaled, real: xl/dstar)
c     t       time (rescaled, real: t/dstar)
c     xs      shift since t zero (0 for spatial cases)
c     nxin    collocation points in x
c     nypin   collocation points in y
c     nzcin   collocation points in z
c     nfzsin  symmetry flag (0: no symmetry)
c     fltype  flow type
c            -3: Asymptotic suction boundary layer
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


      if (scalar.ge.1) then
c
c     Read Prandtl number
c
         m1 = 0.
         pr = 0.
         gr = 0.
         read(12,err=2001) re,pou,xl,zl,t,xs,(pr(i),m1(i),i=1,scalar)

         do i=1,scalar
            if (m1(i).eq.0.or.abs(m1(i)-0.5).lt.1.e-13) then
            else
               if (my_node.eq.0) then
                  write(ioe,*) 'Variation of temperature profile m1'
                  write(ioe,*) 'not implemented.',m1(i),' Scalar no.',i
               end if
               call stopnow(6654)
            end if
         end do
      else
c
c     Ignore Prandtl number and set dummy value
c
         read(12,err=2001) re,pou,xl,zl,t,xs
      end if
      read(12) nxin,nypin,nzcin,nfzsin
      fltype= 0
      rlam  = 0.
      dstar = 0.
      spanv = 0.
      read(12,err=1010) fltype,dstar
 1010 continue
c
c     Write flow type to standard output
c
      if (my_node.eq.0) then
         write(ios,*) 'Reading velocity field'
         write(ios,*) '  filename    : ',trim(namnin)
         write(ios,*) '  m           : ',m
         if (fltype.eq.-3) write(ios,*)
     &        '  fltype      : -3. ',
     &        'Asymptotic suction boundary layer'
         if (fltype.eq.-2) write(ios,*)
     &        '  fltype      : -2. ',
     &        'Temporal Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.-1) write(ios,*)
     &        '  fltype      : -1. Temporal Falkner-Skan boundary layer'
         if (fltype.eq.0) write(ios,*)
     &        '  fltype      : 0. No base flow'
         if (fltype.eq.1) write(ios,*)
     &        '  fltype      : 1. Temporal Poiseuille flow'
         if (fltype.eq.2) write(ios,*)
     &        '  fltype      : 2. Temporal Couette flow'
         if (fltype.eq.3) write(ios,*)
     &        '  fltype      : 3. Temporal Blasius boundary layer'
         if (fltype.eq.4) write(ios,*)
     &        '  fltype      : 4. Spatial Poiseuille flow'
         if (fltype.eq.5) write(ios,*)
     &        '  fltype      : 5. Spatial Couette flow'
         if (fltype.eq.6) write(ios,*)
     &        '  fltype      : 6. Spatial Blasius boundary layer'
         if (fltype.eq.7) write(ios,*)
     &        '  fltype      : 7. Spatial Falkner-Skan boundary layer'
         if (fltype.eq.8) write(ios,*)
     &        '  fltype      : 8. ',
     &        'Spatial Falkner-Skan-Cooke boundary layer'
         if (fltype.eq.9) write(ios,*)
     &        '  fltype      : 9. Spatial parallel boundary layer'
         if (fltype.eq.-20) write(ios,*)
     &        '  fltype      : -20. Temporal buoyancy boundary layer'
         if (fltype.eq.20) write(ios,*)
     &        '  fltype      : 20. Spatial buoyancy boundary layer'
      end if

      if ( fltype.ge.-3.and.fltype.le.-1.or.
     &     fltype.ge.1 .and.fltype.le.9.or.
     &     abs(fltype).eq.20) then
c
c     Correct fltype
c
      else
         write(ioe,*) 'The input file does not contain'
         write(ioe,*) 'the correct type of flow, now:',fltype
         stop
      end if

      if (spat) then
         if ( fltype.eq.-2.or.fltype.eq.-1.or.fltype.eq.1.or.
     &        fltype.eq. 2.or.fltype.eq. 3.or.fltype.eq.-20) then
            write(ioe,*) 'Conflicting variables. Spatial flow but '//
     &           'temporal flow type.'
            write(*,*) 'Change spat in bla.i or use other flow field.'
            stop
         end if
      else
         if ( fltype.eq.4.or.fltype.eq.5.or.fltype.eq.6.or.
     &        fltype.eq.7.or.fltype.eq.8.or.fltype.eq.9.or.
     &        fltype.eq.20) then
            write(ioe,*) 'Conflicting variables. Temporal flow but '//
     &           'spatial flow type.'
            write(*,*) 'Change spat in bla.i or use other flow field.'
            stop
         end if
      end if
c
c     Ensure that dstar is correctly defined in
c     channel and Couette flow cases (=1)
c
      if (fltype.eq.1.or.fltype.eq.2.or.
     &    fltype.eq.4.or.fltype.eq.5) then
          dstar=1.
      end if
c
c     Read additional info for specific flow types
c
      if (fltype.eq.-1) read(12) rlam
      if (fltype.eq.-2) read(12) rlam,spanv
      if (fltype.ge.4.and.fltype.le.6) then
         read(12) bstart,blength
         rlam=0.0
         spanv=0.0
      end if
      if (fltype.ge.7.and.fltype.le.9)
     &     read(12) bstart,blength,rlam,spanv

      if (abs(fltype).eq.20) read(12) (gr(i),i=1,scalar)

c     for ny = 33: 1708 is unstable for k=3.12
c     for Ra=5000 and size 8.66025x15 --> squares

c      gr(1) = 1./8.*pr(1) * 4000.


      if (my_node.eq.0) then
         if (scalar.ge.1) then
            do i=1,scalar
               write(ios,'(4(a,f10.3))')
     &              '   Re          : ',re*dstar,
     &              ' Pr :',pr(i),' m1 :',m1(i),' Gr :',gr(i)
            end do
         else
            write(ios,'(a,f14.5)')
     &           '   Re          : ',re*dstar
         end if
      end if


      nxz=nx/2*mbz
      nxtmp=nxin+nfxd*nxin/2
c
c     Check file info
c
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &     nfzsin.ne.nfzsym) then
         if (my_node.eq.0) then
            write(ioe,*) 'Input file has a size other than program'
            write(ioe,'(a,4i5)') '   File parameters:    ',
     &           nxin,nypin,nzcin,nfzsin
            write(ioe,'(a,4i5)') '   Program parameters: ',
     &           nx,nyp,nzc,nfzsym
            write(ioe,*) '   (nx,nyp,nzc,nfzsym)'
         end if
         if (.not.varsiz) then
            call stopnow(453565)
         end if
         if (nypin.gt.nyp) then
            if (my_node.eq.0) then
               write(*,*) 'Resolution cannot be reduced in y-direction'
            end if
            call stopnow(326576)
         end if
         if (nfzsin.ne.nfzsym) then
            if (my_node.eq.0) then
               write(ioe,*) 'Symmetry cannot be changed'
            end if
            call stopnow(6764)
         end if
      end if
c
c     Prepare FFTs
c
      nxr=min(nx,nxin)
      if (nypin.lt.nyp) then
         call vcosti(nypin,preyin,0)
      end if
      call vcosti(nyp,prey,0)
c
c     Read data from file
c
      do i=1,m
         do zb=1,nzc,mbz
            do z=zb,zb+mbz-1

               if ((nfzsym.eq.0.and.z.gt.(nzcin+1)/2
     &              .and.z.le.nz-nzcin/2)
     &              .or.(nfzsym.eq.1.and.z.gt.nzcin)) then
c
c     If expanding in z-direction and record is new, pad with zeroes
c
                  do y=1,nypin
                     do x=1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                  end do
               else
c
c     Otherwise, read an x-y plane
c
                  do y=1,nypin
                     read(12) (urx(x),x=1,nxr)
                     do x=1,nxr/2
                        boxr(x,z-zb+1,y)=urx(2*x-1)
                        boxi(x,z-zb+1,y)=urx(2*x)
                     end do
c
c     Set boundary velocities to mean in the streamwise direction (0/0-mode)
c     (could be not accurate for varying freestream velocity)
c     It will however be reset in bflow (simulations with fringe), but then
c     according to the value at the inflow (x=0)
c
                     if (z.eq.1.and.y.eq.1    .and.i.eq.1) u0upp=urx(1)
                     if (z.eq.1.and.y.eq.1    .and.i.eq.3) w0upp=urx(1)
                     if (z.eq.1.and.y.eq.nypin.and.i.eq.1) u0low=urx(1)
                     if (z.eq.1.and.y.eq.nypin.and.i.eq.3) w0low=urx(1)

                  end do
c
c     If expanding in x-direction, padding with zeros
c
                  do y=1,nypin
                     do x=nxin/2+1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                  end do
               end if
c
c     If contracting in z-direction then skip records on file
c
               if (z.eq.nz/2.or.nz.eq.1) then
                  do zs=nzc+1,nzcin
                     do y=1,nypin
                        read(12) (urx(x),x=1,nxr)
                     end do
                  end do
               end if

            end do
c
c     If expanding in y-direction then Chebyshev transform, expand and return
c
            if (nypin.lt.nyp) then
               call vchbf(boxr,w3,nypin,nxz,nxz,1,preyin)
               call vchbf(boxi,w3,nypin,nxz,nxz,1,preyin)
               do y=1,nypin
                  do z=1,mbz
                     do x=1,nx/2
                        boxr(x,z,y)=boxr(x,z,y)*(2./real(nypin-1))
                        boxi(x,z,y)=boxi(x,z,y)*(2./real(nypin-1))
                     end do
                  end do
               end do
               do y=nypin+1,nyp
                  do z=1,mbz
                     do x=1,nx/2
                        boxr(x,z,y)=0.0
                        boxi(x,z,y)=0.0
                     end do
                  end do
               end do
               call vchbb(boxr,w3,nyp,nxz,nxz,1,prey)
               call vchbb(boxi,w3,nyp,nxz,nxz,1,prey)
            end if
c
c     This plane belongs to proc: (zb-1)/(nzc/nproc)
c     and is plane number: mod(zb-1,nzc/nproc)+1
c
            node_t = (zb-1)/(nzc/nproc)
            zb_t   = mod(zb-1,nzc/nproc)+1

            if (my_node.eq.node_t) then
c
c     While expanding in z, zero the old oddball modes (zbb)
c
               if (nfzsym.eq.0.and.nzc.gt.nzcin.and.
     &              mod(nzcin,2).eq.0) then
                  zu=nzc-nzcin/2+1
                  zbb=(zu-1)/mbz*mbz+1
                  if (zbb.eq.zb) then
                     do y=1,nyp
                        do x=1,nx/2
                           boxr(x,zu-zb+1,y)=0.0
                           boxi(x,zu-zb+1,y)=0.0
                        end do
                     end do
                  end if
               end if

               if (i.ge.4) then
c
c     Put the scalar into 8
c
                  call putxy(boxr,boxi,zb_t,
     &                 8+pressure+3*(i-4),ur,ui)
               else
c
c     Put the velocities into 1..3
c
                  call putxy(boxr,boxi,zb_t,i,ur,ui)
               end if
            end if

         end do
      end do
c
c     Close file
c
      close(unit=12)
c
c     This is for temporal simulations with derivative bc
c
      du0upp=0.0
c
c     u0upp etc. will be (possibly) overwritten in bflow.f
c
      if (my_node.eq.0) then
         write(ios,*) '  u0low,w0low : ',
     &        u0low,w0low
         write(ios,*) '  u0upp,w0upp : ',
     &        u0upp,w0upp
         write(ios,*) '        spanv : ',spanv
      end if

c
c     change time (or any other parameter) if needed

c      t=0
c      xl = (2*3.1415926/0.303)*dstar
c      xl = 20.736584 *dstar
c      re = 519.15/dstar
c      re = 519.4/dstar
c     critical is 519.4, 0.303
c     unstable at 519

c
c     ----------------------
c      t  =   0.*dstar
c      re = 520./dstar
c     ----------------------
c
c     change domain size
c
      if (1.eq.0) then
c
c     and double the flow in x
c
         xl = 2*xl
         do i=1,memnxyz
            do z=1,memnz
               do y=1,memny
                  do x=2,memnx/2
                     urx(x) = ur(x,y,z,i)
                  end do
                  do x=2,memnx
                     ur(x,y,z,i) = 0.
                  end do
                  do x=2,memnx/2
                     ur(x*2-1,y,z,i) = urx(x)
                  end do
                  
                  do x=2,memnx/2
                     urx(x) = ui(x,y,z,i)
                  end do
                  do x=2,memnx
                     ui(x,y,z,i) = 0.
                  end do
                  do x=2,memnx/2
                     ui(x*2-1,y,z,i) = urx(x)
                  end do
               end do
            end do
         end do
      end if
      if (1.eq.0) then
c
c     and double the flow in z
c
         zl = 2*zl
         do i=1,memnxyz
            do y=1,memny
               do x=1,memnx
                  do z=2,memnz
                     urz(z) = ur(x,y,z,i)
                  end do
                  do z=2,memnz
                     ur(x,y,z,i) = 0.
                  end do
                  do z=2,memnz/4
                     ur(x,y,z*2-1,i) = urz(z)
                  end do
                  do z=1,memnz/4
                     ur(x,y,memnz-memnz/2+z*2-1,i) = 
     &                    urz(memnz-memnz/4+z)
                  end do
                  ur(x,y,memnz/2+1,i) = 0.
                  
                  do z=2,memnz
                     urz(z) = ui(x,y,z,i)
                  end do
                  do z=2,memnz
                     ui(x,y,z,i) = 0.
                  end do
                  do z=2,memnz/4
                     ui(x,y,z*2-1,i) = urz(z)
                  end do
                  do z=1,memnz/4
                     ui(x,y,memnz-memnz/2+z*2-1,i) = 
     &                    urz(memnz-memnz/4+z)
                  end do
                  ui(x,y,memnz/2+1,i) = 0.
               end do
            end do
         end do
      end if
c
c     impose symmetry
c
      if (1.eq.0) then
         write(*,*) 'impose spanwise symmetry wrt z=0'
         do y=1,memny
            do z=2,memnz/2
               do x=1,memnx
                  ur(x,y,memnz+2-z,1) = ur(x,y,z,1)
                  ui(x,y,memnz+2-z,1) = ui(x,y,z,1)
                  ur(x,y,memnz+2-z,2) = ur(x,y,z,2)
                  ui(x,y,memnz+2-z,2) = ui(x,y,z,2)
                  ur(x,y,memnz+2-z,3) = -ur(x,y,z,3)
                  ui(x,y,memnz+2-z,3) = -ui(x,y,z,3)
               end do
            end do
         end do
      end if

      return

 2001 continue
      if (my_node.eq.0) then
         write(ioe,*) 'Error reading file header, line 1:'
         write(ioe,*) 're = ',re
         write(ioe,*) 'pou= ',pou
         write(ioe,*) 'xl = ',xl
         write(ioe,*) 'zl = ',zl
         write(ioe,*) 't  = ',t
         write(ioe,*) 'xs = ',xs
         do i=1,scalar
            write(ioe,*) 'pr = ',pr(i),i
            write(ioe,*) 'm1 = ',m1(i),i
         end do
      end if
      call stopnow(2001)

      end subroutine rdiscbl



      subroutine swap(a,b)
c
c     Swaps the content of two reals
c
      implicit none
      real a,b,c
      c=b
      b=a
      a=c
      end subroutine swap



      subroutine mirror_u(ur,ui)
c
c     Mirrors the velocity field ur,ui in various directions.
c     This routine acts on the velocities ur/ui(..,..,..,1-3) only.
c     Note that currently only the complete mirroring is tested!
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer x,y,z

      if (1.eq.1) then
         write(*,*) 'impose spanwise mirror wrt z=0'
         do y=1,memny
            do z=2,memnz/2
               do x=1,memnx
                  call swap(ur(x,y,memnz+2-z,1),ur(x,y,z,1))
                  call swap(ui(x,y,memnz+2-z,1),ui(x,y,z,1))
                  call swap(ur(x,y,memnz+2-z,2),ur(x,y,z,2))
                  call swap(ui(x,y,memnz+2-z,2),ui(x,y,z,2))
                  call swap(ur(x,y,memnz+2-z,3),ur(x,y,z,3))
                  call swap(ui(x,y,memnz+2-z,3),ui(x,y,z,3))
               end do
            end do
         end do
         ur(:,:,:,3) = -ur(:,:,:,3)
         ui(:,:,:,3) = -ui(:,:,:,3)
      end if
      if (1.eq.1) then
         write(*,*) 'impose streamwise mirror wrt x=xlen/2'
         do y=1,memny
            do z=2,memnz/2
               do x=1,memnx
                  call swap(ur(x,y,memnz+2-z,1),ur(x,y,z,1))
                  call swap(ui(x,y,memnz+2-z,1),ui(x,y,z,1))
                  call swap(ur(x,y,memnz+2-z,2),ur(x,y,z,2))
                  call swap(ui(x,y,memnz+2-z,2),ui(x,y,z,2))
                  call swap(ur(x,y,memnz+2-z,3),ur(x,y,z,3))
                  call swap(ui(x,y,memnz+2-z,3),ui(x,y,z,3))
               end do
            end do
            do z=1,memnz
               do x=1,memnx
                  ur(x,y,z,1)=-ur(x,y,z,1)
                  ui(x,y,z,2)=-ui(x,y,z,2)
                  ui(x,y,z,3)=-ui(x,y,z,3)
               end do
            end do
         end do
      end if
      if (1.eq.1) then
         write(*,*) 'impose vertical mirror wrt y=0'
         do y=1,(memny-1)/2
            do z=1,memnz
               do x=1,memnx
                  call swap(ur(x,memny+1-y,z,1),ur(x,y,z,1))
                  call swap(ui(x,memny+1-y,z,1),ui(x,y,z,1))
                  call swap(ur(x,memny+1-y,z,2),ur(x,y,z,2))
                  call swap(ui(x,memny+1-y,z,2),ui(x,y,z,2))
                  call swap(ur(x,memny+1-y,z,3),ur(x,y,z,3))
                  call swap(ui(x,memny+1-y,z,3),ui(x,y,z,3))
               end do
            end do
         end do
         do y=1,memny
            do z=1,memnz
               do x=1,memnx
                  ur(x,y,z,2) = -ur(x,y,z,2)
                  ui(x,y,z,2) = -ui(x,y,z,2)
               end do
            end do
         end do
      end if

      end subroutine mirror_u


      subroutine mirror_v(ur,ui,puw)
c
c     Mirrors the velocity field ur,ui in various directions.
c     This routine acts on the vorticities and prhs ur/ui(..,..,..,4-7).
c     Also, the mean puw is adjusted.
c     Note that currently only the complete mirroring is implemented!
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      integer x,y,z
      real puw(ny,2+scalar)
c
c     in z
c
      do y=1,memny
         do z=2,memnz/2
            do x=1,memnx
               call swap(ur(x,y,memnz+2-z,4),ur(x,y,z,4))
               call swap(ui(x,y,memnz+2-z,4),ui(x,y,z,4))
               call swap(ur(x,y,memnz+2-z,5),ur(x,y,z,5))
               call swap(ui(x,y,memnz+2-z,5),ui(x,y,z,5))
            end do
         end do
      end do
c
c     in x
c
      do y=1,memny
         do z=2,memnz/2
            do x=1,memnx
               call swap(ur(x,y,memnz+2-z,4),ur(x,y,z,4))
               call swap(ui(x,y,memnz+2-z,4),ui(x,y,z,4))
               call swap(ur(x,y,memnz+2-z,5),ur(x,y,z,5))
               call swap(ui(x,y,memnz+2-z,5),ui(x,y,z,5))
            end do
         end do
         do z=1,memnz
            do x=1,memnx
               ur(x,y,z,4)=-ur(x,y,z,4)
               ur(x,y,z,5)=-ur(x,y,z,5)
            end do
         end do
      end do


      do y=1,memny,2
         do z=1,memnz
            do x=1,memnx
               ur(x,y,z,6)=-ur(x,y,z,6)
               ui(x,y,z,7)=-ui(x,y,z,7)
            end do
         end do
      end do
      do y=2,memny,2
         do z=1,memnz
            do x=1,memnx
               ui(x,y,z,6)=-ui(x,y,z,6)
               ur(x,y,z,7)=-ur(x,y,z,7)
            end do
         end do
      end do
c
c     in y
c      
      do y=1,(memny-1)/2
         do z=1,memnz
            do x=1,memnx
               call swap(ur(x,memny+1-y,z,4),ur(x,y,z,4))
               call swap(ui(x,memny+1-y,z,4),ui(x,y,z,4))
               call swap(ur(x,memny+1-y,z,5),ur(x,y,z,5))
               call swap(ui(x,memny+1-y,z,5),ui(x,y,z,5))
            end do
         end do
      end do

      ur(:,:,:,4) = -ur(:,:,:,4)
      ui(:,:,:,4) = -ui(:,:,:,4)
      ur(:,:,:,5) = -ur(:,:,:,5)
      ui(:,:,:,5) = -ui(:,:,:,5)
c
c     Now also for the mean value
c
      do y=1,ny,2
         puw(y,1) = 0.
         puw(y,2) = 0.
      end do

      end subroutine mirror_v


      subroutine symmetrise_u(ur,ui)
c
c     Symmetrises the field ur,ui using the symmetries in mirror.
c     Here, components 1-3 (velocities) are done.
c     Note that the work array vr,vi is only defined in this subroutine.
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real vr(memnx,memny,memnz,memnxyz),vi(memnx,memny,memnz,memnxyz)
      integer i

      vr = ur
      vi = ui
      call mirror_u(vr,vi)
      ur = 0.5*(ur+vr)
      ui = 0.5*(ui+vi)
    
      end subroutine symmetrise_u


      subroutine symmetrise_v(ur,ui,puw)
c
c     Symmetrises the field ur,ui using the symmetries in mirror.
c     Here, components 4,5 (vorticities) and 6,7 (partial RHS) are done.
c     Note that the work array vr,vi is only defined in this subroutine.
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real vr(memnx,memny,memnz,memnxyz),vi(memnx,memny,memnz,memnxyz)
      real puw(ny,2+scalar)

      vr = ur
      vi = ui
      call mirror_v(vr,vi,puw)
      ur = 0.5*(ur+vr)
      ui = 0.5*(ui+vi)

      end subroutine symmetrise_v
