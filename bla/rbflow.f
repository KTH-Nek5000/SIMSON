c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rbflow(bu1,bu2,nambfl,u0low,w0low,
     &     prex,wbr,wbi,boxr,boxi,w3)
c
c     Reads the spanwise-averaged flow from a standard 3D flow file
c
      implicit none

      include 'par.f'

      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      character*80 nambfl
      real u0low,w0low
      real prex(nxp+15)
      real wbr((nxp/2+1),nyp),wbi((nxp/2+1),nyp)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real w3(nx/2,mbz,nyp)
      real gr(scalar)

      integer fltype
      real u0lowb,w0lowb
      real re,xl,zl,t,xs,dstar,bstart,blength,spanv,rlam
      logical pou
      real urx(nx)

      integer x,y,z,i,zb,zs,nxz,nxr
      integer nxin,nypin,nzcin,nfzsin
      real preyin(2*nyp+15),prey(2*nyp+15)

      open(unit=112,file=nambfl,status='old',form='unformatted')
c
c     Read file header
c
      read(112) re,pou,xl,zl,t,xs
      read(112) nxin,nypin,nzcin,nfzsin
      read(112) fltype,dstar

      if (fltype.ne.3.and.(fltype.lt.4.or.fltype.gt.8)
     &     .and.fltype.ne.20) then
         write(ioe,*) 'The base flow file does not contain'
         write(ioe,*) 'the correct type of flow'
         stop
      end if

      bstart=0.0
      blength=0.0
      rlam=0.0
      spanv=0.0
      if (fltype.ge.4.and.fltype.le.6) read(112) bstart,blength
      if (fltype.eq.-1) read(112) rlam
      if (fltype.eq.-2) read(112) rlam,spanv
      if (fltype.ge.7.and.fltype.le.9)
     &     read(112) bstart,blength,rlam,spanv

      if (fltype.eq.20) read(112) (gr(i),i=1,scalar)

      nxz=nx/2*mbz

c      px=0.
c      if (pou) px=-2./re

      if (nypin.gt.nyp) then
         write(*,*) 'The base flow cannot have a higher y-resolution'
         stop
      end if
      nxr=min(nx,nxin)
      if (nypin.lt.nyp) then
         call vcosti(nypin,preyin,0)
      end if
      call vcosti(nyp,prey,0)

      do i=1,3+scalar
         do zb=1,nzc,mbz
            do z=zb,zb+mbz-1
c
c     If expanding in z-direction and record is new pad with zeroes
c
               if ((nfzsym.eq.0.and.z.gt.(nzcin+1)/2
     &              .and.z.le.nz-nzcin/2)
     &              .or.(nfzsym.eq.1.and.z.gt.nzcin)) then
                  do y=1,nypin
                     do x=1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                  end do
               else
                  do y=1,nypin
                     read(112) (urx(x),x=1,nxr)
                     do x=1,nxr/2
                        boxr(x,z-zb+1,y)=urx(2*x-1)
                        boxi(x,z-zb+1,y)=urx(2*x)
                     end do
                     do x=nxin/2+1,nx/2
                        boxr(x,z-zb+1,y)=0.0
                        boxi(x,z-zb+1,y)=0.0
                     end do
                     if (z.eq.1.and.y.eq.nypin.and.i.eq.1) u0lowb=urx(1)
                     if (z.eq.1.and.y.eq.nypin.and.i.eq.3) w0lowb=urx(1)
                  end do
               end if
c
c     If contracting in z-direction then skip records on file
c
               if (z.eq.nz/2.or.nz.eq.1) then
                  do zs=nzc+1,nzcin
                     do y=1,nypin
                        read(112) (urx(x),x=1,nxr)
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
            if (zb.eq.1) then
               do y=1,nyp
                  do x=1,nx/2
                     bu1(x,y,i)=boxr(x,1,y)
                     bu2(x,y,i)=boxi(x,1,y)
                  end do
               end do
            end if
         end do
      end do
      close(unit=112)
c
c     Adjust the lower wall velocity to that of the starting field
c
      do y=1,nyp
         bu1(1,y,1)=bu1(1,y,1)+u0low-u0lowb
         bu1(1,y,3)=bu1(1,y,3)+w0low-w0lowb
      end do
c
c     Transform to physical space in x-direction
c
      do i=1,3+scalar
         do y=1,nyp
            do x=nx/2+1,nxp/2+1
               bu1(x,y,i)=0.0
               bu2(x,y,i)=0.0
            end do
         end do
         call vrfftb(bu1(1,1,i),bu2(1,1,i),wbr,wbi,nxp,nyp,
     &        1,nxp/2+1,prex)
      end do

      end subroutine rbflow
