c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wplbl_serial(npl,nampl,tpl,cpl,it,b2r,b2i,uxz,uxy,
     &     uyz,wr,wi,ur,ui,mpl,re,xl,zl,t,xs,dstar,fltype,
     &     prexn,prezn,presn,prean,u0low,w0low,alfa,zs,beta,wplblodd)
c
c     Writes npl planes to the same number of files
c     note that for MPI data transpose is necessary.
c
      implicit none

      include 'par.f'

      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real b2r(nxp/2+1,mby,nzd),b2i(nxp/2+1,mby,nzd)
      real wr(nxp/2+1,mby,nzd),wi(nxp/2+1,mby,nzd)
      integer nzn
      parameter (nzn=nz-nfzsym*(nz/2-1))
      real uxz(nx,nzn),uxy(nx,nyp),uyz(nyp,nzn),u0low,w0low
      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real alfa(nx/2*mbz),zs,beta(nz)

      integer npl,mpl,tpl(mpl,3)
      real cpl(mpl)
      character*80 nampl(mpl)
      integer it,fltype
      real t,re,xl,zl,xs,dstar
      integer wplblodd(mpl)
      real pi
      parameter (pi = 3.1415926535897932385)
      integer i,x,y,z,ivar,yp,yb,zp,xp,j
      real cx(2,mpl),cxx(2),cy(2,mpl),cyy(2),cz(mpl),czz
      integer cxi(2),cyi(2),czi,co
      real data(nx*ny)
      logical fapp
      integer,save :: firsttime=0

      write(ios,*) '** write planes (serial) **',t/dstar

      do i=1,npl
c
c     Variable to save
c
         ivar=tpl(i,2)
c
c     If first iteration, open files and write header
c     If file exists, do not overwrite it, append the data
c
c         if (it.eq.1) then
         if (firsttime.eq.0) then
            inquire(file=nampl(i),exist=fapp)
            if(fapp) then
               write(ios,*) '*** planes appended to file'
               open(unit=30+i,file=nampl(i),form='unformatted',
     &              position='append')
            else
               open(unit=30+i,file=nampl(i),form='unformatted')
            end if
c
c     Write headers
c
            if (tpl(i,1).eq.1) then
c
c     xy plane
c
               cpl(i)=cpl(i)*dstar
               tpl(i,3)=int(cpl(i)/zl*real(nz)+.5)+nz/2+1
               cpl(i)=zl/real(nz)*real(tpl(i,3)-nz/2-1)/dstar
               if (nfzsym.eq.1.and.cpl(i).gt.0.) tpl(i,3)=nz+2-tpl(i,3)
            else if (tpl(i,1).eq.2) then
c
c     xz plane
c
c     The formula for cpl in the Couette and Channel geometries
c     allows for negative y-coordinates

               if (fltype.eq.1.or.fltype.eq.2.or.
     &              fltype.eq.4.or.fltype.eq.5) then

                  cpl(i)=cpl(i)*dstar
                  tpl(i,3)=int(acos(cpl(i))/pi*real(nyp-1)+1.5)
                  cpl(i)=(cos(pi*real(tpl(i,3)-1)/real(nyp-1)))/dstar

               else

                  cpl(i)=cpl(i)*dstar-1.
                  tpl(i,3)=int(acos(cpl(i))/pi*real(nyp-1)+1.5)
                  cpl(i)=(cos(pi*real(tpl(i,3)-1)/real(nyp-1))+1.)/dstar

               end if

            else if (tpl(i,1).eq.3) then
c
c     yz plane
c
               cpl(i)=cpl(i)*dstar
               tpl(i,3)=int(cpl(i)/xl*real(nx)+.5)+nx/2+1
               tpl(i,3)=tpl(i,3)-(tpl(i,3)-1)/nx*nx
               wplblodd(i)=mod(tpl(i,3),2)
               tpl(i,3)=(tpl(i,3)+1)/2
               cpl(i)=int(cpl(i)/xl*real(nx)+.5)*xl/real(nx)
            else
               call stopnow(454532)
            end if
c
c     The data will be shifted so that xs,zs=0
c
            if(fapp.eqv..false.) then
               write(30+i) re,.false.,xl,zl,t,0.
               write(30+i) nx,nyp,nzc,nfzsym
               write(30+i) tpl(i,1),ivar,cpl(i),fltype,dstar
            end if
         end if
c
c     xy planes
c
         if (tpl(i,1).eq.1) then
c
c     Note that subplane saving is implemented (only for xy planes)
c     this way, the plane file has two additional records at plane no. 1
c
            zp=tpl(i,3)
            do yb=1,nyp,mby
               call getxz(b2r,b2i,yb,ivar,0,ur,ui)
c
c     No shift for Couette flow
c
               if (fltype.ne.2.and.fltype.ne.5) then
                  call xzsh(b2r,b2i,xs,zs,alfa,beta,yb)
                  if (ivar.eq.1) then
                     do y=1,min(mby,nyp-yb+1)
                        b2r(1,y,1)=b2r(1,y,1)-u0low
                     end do
                  end if
                  if (ivar.eq.3) then
                     do y=1,min(mby,nyp-yb+1)
                        b2r(1,y,1)=b2r(1,y,1)-w0low
                     end do
                  end if
               end if


               do y=1,min(mby,nyp-yb+1)
                  if (nfzsym.eq.0) then
                     call vcfftb(b2r(1,y,1),b2i(1,y,1),wr,wi,nz,
     &                    nx/2,(nxp/2+1)*mby,1,prezn)
                  else
                     if (ivar.le.2) then
                        call vcffts(b2r(1,y,1),b2i(1,y,1),wr,
     &                       nz/2+1,nx/2,(nxp/2+1)*mby,1,presn)
                     else
                        call vcftab(b2r(1,y,2),b2i(1,y,2),
     &                       wr,wi,nz/2-1,nx/2,(nxp/2+1)*mby,1,prean)
                     end if
                  end if
                  call vrfftb(b2r(1,y,zp),b2i(1,y,zp),wr,wi,
     &                 nx,1,1,nxp/2+1,prexn)
               end do
               do y=1,min(mby,nyp-yb+1)
                  do x=1,nx/2
                     uxy(2*x-1,yb+y-1)=b2r(x,y,zp)
                     uxy(2*x,yb+y-1)=b2i(x,y,zp)
                  end do
               end do
            end do

            if (0.eq.1) then
c
c     Use subplane saving
c     coordinates for subplane 1
c
               cx(1,1)=250.
               cx(2,1)=600.
               cy(1,1)=0.
               cy(2,1)=10.
               cz(1)=10.
c
c     Coordinates for subplane 2
c
               cx(1,2)=250.
               cx(2,2)=600.
               cy(1,2)=0.
               cy(2,2)=10.
               cz(2)=7.
c
c     Coordinates for subplane 3
c
               cx(1,3)=250.
               cx(2,3)=600.
               cy(1,3)=0.
               cy(2,3)=10.
               cz(3)=7.

               do x=1,2
c
c     x-coordinate
c
                  cxi(x)=int(cx(x,i)*dstar/xl*real(nx)+.5)+nx/2+1
                  cxi(x)=cxi(x)-(cxi(x)-1)/nx*nx
                  cxx(x)=int(cx(x,i)*dstar/xl*real(nx)+.5)*
     &                 xl/dstar/real(nx)
c                write(ios,*) 'x ',cx(x,i),cxx(x),cxi(x)
c
c     y-coordinate
c
                  cyy(x)=cy(x,i)*dstar-1.
                  cyi(x)=int(acos(cyy(x))/pi*real(nyp-1)+1.5)
                  cyy(x)=(cos(pi*real(cyi(x)-1)/real(nyp-1))+1.)/dstar
c                write(ios,*) cy(x),cyy(x),cyi(x)
               end do
c
c     z-coordinate
c
               czi=int(cz(i)*dstar/zl*real(nz)+.5)+nz/2+1
               czz=zl/real(nz)*real(czi-nz/2-1)/dstar
c             write(ios,*) cz,czz,czi
               if (cyi(2).lt.cyi(1)) then
c
c     Swap y coordinates
c
                  x=cyi(1)
                  cyi(1)=cyi(2)
                  cyi(2)=x
               end if

               if (it.eq.1) then
c
c     In the first iteration, write out coordinates
c
                  do x=1,2
                     write(30+i) cxx(x),cxi(x),cyy(x),cyi(x)
                  end do
                  write(30+i) czz,czi
               end if

               co=0
               if (cxi(1).lt.cxi(2)) then
c
c     Normal order...
c
                  do x=cxi(1),cxi(2)
                     do y=cyi(1),cyi(2)
                        co=co+1
                        data(co)=uxy(x,y)
                     end do
                  end do
               else
c
c     Wrapped over edge
c
                  do x=cxi(1),nx
                     do y=cyi(1),cyi(2)
                        co=co+1
                        data(co)=uxy(x,y)
                     end do
                  end do
                  do x=1,cxi(2)
                     do y=cyi(1),cyi(2)
                        co=co+1
                        data(co)=uxy(x,y)
                     end do
                  end do
               end if
               if (it.eq.1) then
c
c     In the first iteration, write out number of data
c
                  write(30+i) co
               end if
c
c     And write out data
c
               write(30+i) t
               write(30+i) (data(j),j=1,co)
            else
c
c     Write out data (whole plane)
c
               write(30+i) t,0.
               write(30+i) uxy
            end if
         end if
c
c     xz planes
c
         if (tpl(i,1).eq.2) then
            write(30+i) t,0.
            yp=tpl(i,3)
            yb=(yp-1)/mby*mby+1
            call getxz(b2r,b2i,yb,ivar,0,ur,ui)
c
c     No shift for Couette flow
c
               if (fltype.ne.2.and.fltype.ne.5) then
                  call xzsh(b2r,b2i,xs,zs,alfa,beta,yb)
                  if (ivar.eq.1) b2r(1,yp-yb+1,1)=b2r(1,yp-yb+1,1)-u0low
                  if (ivar.eq.3) b2r(1,yp-yb+1,1)=b2r(1,yp-yb+1,1)-w0low
               end if
            if (nfzsym.eq.0) then
               call vcfftb(b2r(1,yp-yb+1,1),b2i(1,yp-yb+1,1),wr,wi,
     &              nz,nx/2,(nxp/2+1)*mby,1,prezn)
            else
               if (ivar.le.2) then
                  call vcffts(b2r(1,yp-yb+1,1),b2i(1,yp-yb+1,1),wr,
     &                 nz/2+1,nx/2,(nxp/2+1)*mby,1,presn)
               else
                  call vcftab(b2r(1,yp-yb+1,2),b2i(1,yp-yb+1,2),wr,wi,
     &                 nz/2-1,nx/2,(nxp/2+1)*mby,1,prean)
               end if
            end if
            call vrfftb(b2r(1,yp-yb+1,1),b2i(1,yp-yb+1,1),wr,wi,
     &           nx,nzpc,1,mby*(nxp/2+1),prexn)
            do z=1,nzn
               do x=1,nx/2
                  uxz(2*x-1,z)=b2r(x,yp-yb+1,z)
                  uxz(2*x,z)=b2i(x,yp-yb+1,z)
               end do
            end do
            write(30+i) uxz
         end if
c
c     yz planes
c
         if (tpl(i,1).eq.3) then
c           write(ios,*) wplblodd(i)
            write(30+i) t,0.
            xp=tpl(i,3)
            do yb=1,nyp,mby
               call getxz(b2r,b2i,yb,ivar,0,ur,ui)
c
c     No shift for Couette flow
c
               if (fltype.ne.2.and.fltype.ne.5) then
                  call xzsh(b2r,b2i,xs,zs,alfa,beta,yb)
                  if (ivar.eq.1) then
                     do y=1,min(mby,nyp-yb+1)
                        b2r(1,y,1)=b2r(1,y,1)-u0low
                     end do
                  end if
                  if (ivar.eq.3) then
                     do y=1,min(mby,nyp-yb+1)
                        b2r(1,y,1)=b2r(1,y,1)-w0low
                     end do
                  end if
               end if
               do y=1,min(mby,nyp-yb+1)
                  if (nfzsym.eq.0) then
                     call vcfftb(b2r(1,y,1),b2i(1,y,1),wr,wi,nz,
     &                    nx/2,(nxp/2+1)*mby,1,prezn)
                  else
                     if (ivar.le.2) then
                        call vcffts(b2r(1,y,1),b2i(1,y,1),wr,
     &                       nz/2+1,nx/2,(nxp/2+1)*mby,1,presn)
                     else
                        call vcftab(b2r(1,y,2),b2i(1,y,2),
     &                       wr,wi,nz/2-1,nx/2,(nxp/2+1)*mby,1,prean)
                     end if
                  end if
                  call vrfftb(b2r(1,y,1),b2i(1,y,1),wr,wi,
     &                 nx,nzpc,1,(nxp/2+1)*mby,prexn)
               end do
               do y=1,min(mby,nyp-yb+1)
                  do z=1,nzn
                     uyz(yb+y-1,z)=
     &                    b2r(xp,y,z)*wplblodd(i)-
     &                    b2i(xp,y,z)*(wplblodd(i)-1)
                  end do
               end do
            end do
            write(30+i) uyz
         end if
      end do

      firsttime = 1


      end subroutine wplbl_serial
