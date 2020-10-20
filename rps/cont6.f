c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine cont6(nx,ny,u1,xy,mxs,mxe,mys,mye,mcont,slev,elev,
     &     xaxis,yaxis,heada,headb,mbox,jhard,lraw,gl,gu,ipln,ifil,
     &     umin,umax,nothing)
c
c     Calculates contours of values in an array u1 with associated
c     coordinates xy
c     plotting is done from mxs to mxe and from mys to mye
c     make mcont contours from slev to elev if mcont positive,
c     if mcont is negative let contxx choose its own levels
c     max abs(mcont)
c
c     Originally by Dan Henningson
c     revised by Anders Lundbladh 890901
c     revised by Henningson spring 91
c     revised by Lundbladh 910515
c     and again 910605
c     changed to cont6 920909
c
      implicit none

      integer naddim,nxdim,ihard,nx,ny,ic,is,isnext,idash,nd
      integer i,mcont,ncont,ncopt,mbox,ifil
      real xmin,xmax,ymin,ymax
      real elev,slev
      parameter(naddim=50000,nxdim=4*naddim)
      real u1(nx,ny),xy(nx,ny,2)
      real acont(naddim),xcont(nxdim),ycont(nxdim)
      integer nlev(naddim),work(naddim),nad(naddim)
      integer jhard,ipln,mxs,mxe,mys,mye
      character text*80
      character*(*) xaxis,yaxis,heada,headb
      real gl,gu
      logical lraw
      common /hard/ihard
      real umin,umax
      logical nothing

      logical wshow

      ihard=jhard

      if (lraw) then
         wshow=ipln.eq.1
         call raw(u1,xy,mxs,mxe,mys,mye,nx,ny,mbox,gl,gu,wshow,ifil,
     &        umin,umax,nothing)
      else
         if (nothing) then
         else
c
c     Contour plots
c
            if (mcont.gt.1) then
               do i=1,mcont
                  acont(i)=slev+(elev-slev)*real(i-1)/real(mcont-1)
               end do
            end if
            if (mcont.eq.1) then
               acont(1)=slev
            end if
            ncopt=2
            if (mcont.lt.0) ncopt=1
            ncont=abs(mcont)
c
            call contxx(nx,ny,mxs,mxe,mys,mye,xy,u1,ncopt,ncont,acont,
     *           nxdim,xcont,ycont,naddim,nad,nlev,naddim,work)
c
            write(text,1222) ncont,acont(1),acont(ncont),acont(2)-
     &           acont(1)
 1222       format(i3,' levels, start:',f11.8,' end:',f10.8,
     &           ' spacing:',f10.8,'$')
            if (abs(acont(1)).GT..1.or. abs(acont(ncont)).GT..1)
     &           write(text,1223) ncont,acont(1),acont(ncont),acont(2)-
     &           acont(1)
 1223       format(i3,' levels, start:',f11.4,' end:',f10.4,
     &           ' spacing:',f10.4,'$')
            xmin=xy(mxs,mys,1)
            xmax=xy(mxe,mye,1)
            ymin=xy(mxs,mys,2)
            ymax=xy(mxe,mye,2)
            call xterm1
            call plarea(xmin,xmax,ymin,ymax,mbox)
            call axis(xmin,0.,xmax,ymin,0.,ymax)
            call head1(heada)
            call head2(headb)
            call head3(text)
            call taxis(xaxis,yaxis)
c
            isnext=1
            nd=nad(1)
            do ic=2,nd+1
               is=isnext
               isnext=nad(ic)
               idash=1
               if (acont(nlev(ic-1)).lt.0) idash=2
               call curve(xcont(is),ycont(is),isnext-is,idash)
            end do
            call plend
         end if
      end if

      return

      end subroutine cont6
