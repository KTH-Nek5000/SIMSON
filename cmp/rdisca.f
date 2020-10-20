c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rdisca(ur,ui,u0low,u0upp,w0low,w0upp,re,mflux,
     &     xl,zl,t,xs,dstar,fltype,bstart,blength,rlam,spanv,
     &     namnin,varsiz,m,ifile,boxr,boxi,w3,urx)
c
c     Reads m variables from file namnin and puts into ur as field ifile
c     reads both channel flow and boundary layer flow
c
      implicit none

      include 'par.f'

      character*32 namnin
      integer m,fltype,ifile
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real u0low,u0upp,w0low,w0upp
      real re,mflux,xl,zl,t,xs,dstar
      real bstart,blength,rlam,spanv
      logical pou,varsiz
      real urx(nx)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real w3(nx/2,mbz,nyp)

      integer x,y,z,i,zb,zs,nxz,nxr,zu,ic
      integer nxin,nypin,nzcin,nfzsin
      real preyin(2*nyp+15),prey(2*nyp+15)
c
c     Open input file, read file header
c
      open(unit=12,file=namnin,status='old',form='unformatted')
      read(12) re,pou,xl,zl,t,xs
      read(12) nxin,nypin,nzcin,nfzsin
      fltype=0
      dstar=0.0
      read(12,err=1010) fltype,dstar
 1010 if (fltype.eq.1) dstar=1.
      if (fltype.eq.2) dstar=1.
      if (fltype.eq.4) dstar=1.
      if (fltype.eq.5) dstar=1.
      if (fltype.eq.-1) read(12) rlam
      if (fltype.eq.-2) read(12) rlam,spanv
      if (fltype.eq.6) then
         read(12) bstart,blength
         rlam = 0.
         spanv =0.
        end if
      if (fltype.ge.7)
     &     read(12) bstart,blength,rlam,spanv
      nxz=nx/2*mbz
c
c     Check file info
c
      if (nxin.ne.nx.or.nypin.ne.nyp.or.nzcin.ne.nzc.or.
     &    nfzsin.ne.nfzsym) then
         write(*,*) 'Input file has a size other than program'
         write(*,*) 'File parameters nx,nyp,nzc,nfzsym',
     &        nxin,nypin,nzcin,nfzsin
         write(*,*) 'Program parameters               ',
     &        nx,nyp,nzc,nfzsym
         if (.not.varsiz) stop
         if (nypin.gt.nyp) then
            write(*,*) 'Resolution cannot be reduced in y-direction'
            stop
         end if
         if (nfzsin.ne.nfzsym) then
            write(*,*) 'Symmetry cannot be changed'
            stop
         end if
      end if
      nxr=min(nx,nxin)
      call vcosti(nypin,preyin,0)
      call vcosti(nyp,prey,0)
c
      do i=1+(ifile-1)*m,ifile*m
         do zb=1,nzc,mbz
            do z=zb,zb+mbz-1
               ic=z+(i-1-(ifile-1)*m)*nzc
               if (ic.eq.nzc*m/100+1) write(*,*) 'Read 1%'
               if (ic.eq.nzc*m/10+1) write(*,*) 'Read 10%'
               if (ic.eq.nzc*m/2+1) write(*,*) 'Read 50%'
               if (ic.eq.nzc*m*9/10+1) write(*,*) 'Read 90%'
c
c     If expanding in z-direction and record is new pad with zeroes
c
               if ((nfzsym.eq.0.and.
     &              z.gt.(nzcin+1)/2.and.z.le.nz-nzcin/2)
     &              .or.(nfzsym.eq.1.and.z.gt.nzcin)) then
                  do y=1,nypin
                     do x=1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                  end do
               else
c
c     Ok here we read an x-y plane
c
                  do y=1,nypin
                     read(12) (urx(x),x=1,nxr)
                     do x=1,nxr/2
                        boxr(x,z-zb+1,y)=urx(2*x-1)
                        boxi(x,z-zb+1,y)=urx(2*x)
                     end do
                  end do
                  do y=1,nypin
                     do x=nxin/2+1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                  end do
                  if (z.eq.1.and.i.eq.1) u0upp=boxr(1,1,1)
                  if (z.eq.1.and.i.eq.3) w0upp=boxr(1,1,1)
                  if (z.eq.1.and.i.eq.1) u0low=boxr(1,1,nypin)
                  if (z.eq.1.and.i.eq.3) w0low=boxr(1,1,nypin)
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
c     If expanding in y-direction then chebyshev transform, expand and return
c
            if (nypin.lt.nyp) then
               call vchbf(boxr,w3,nypin,nxz,nxz,1,preyin)
               call vchbf(boxi,w3,nypin,nxz,nxz,1,preyin)
               do y=1,nypin
                  do z=1,mbz
                     do x=1,nx/2
                        boxr(x,z,y)=boxr(x,z,y)*(2./float(nypin-1))
                        boxi(x,z,y)=boxi(x,z,y)*(2./float(nypin-1))
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
            call putxy(boxr,boxi,zb,i,ur,ui)
            if (zb.eq.1.and.i.eq.1) then
               call vchbf(boxr,w3,nyp,nxz,nxz,1,prey)
               mflux=0.
               do y=1,ny,2
                  mflux=mflux-boxr(1,1,y)*2./float((y-1)**2-1)*
     &                 (2./float(nyp-1))
               end do
               mflux=mflux*zl
            end if
         end do
c
c     Zero old odd ball
c
         if (nfzsym.eq.0.and.nzc.gt.nzcin.and.mod(nzcin,2).eq.0) then
            zu=nzc-nzcin/2+1
            zb=(zu-1)/mbz*mbz+1
            call getxy(boxr,boxi,zb,i,ur,ui)
            do y=1,nyp
               do x=1,nx/2
                  boxr(x,zu-zb+1,y)=0.0
                  boxi(x,zu-zb+1,y)=0.0
               end do
            end do
            call putxy(boxr,boxi,zb,i,ur,ui)
         end if
      end do
      close(unit=12)

      return

      end subroutine rdisca
