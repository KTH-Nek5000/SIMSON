c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine lenbox(ur,ui,bu1,bu2,bucorr,xs,x0,
     &     fbla,rlam,spanv,re,nxtmp,xl,xlold,
     &     dstar,bstart,blength,tabfre,namfre,u0low,w0low,
     &     fltype,fend,fstart,prex,prextm,prez,pres,prea,prey,wbr,wbi,
     &     pxz,w,my_node,eta,dybla,nbla,m1,suction,asbl,vsuc)
c
c     Lengthens the computational box in the streamwise direction by
c     adding a region with only baseflow before start of blending or
c     fringe region
c
      implicit none

      include 'par.f'

      logical tabfre,sym
      character*80 namfre
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
      real bu1(nxp/2+1,nyp,3+scalar),bu2(nxp/2+1,nyp,3+scalar)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real prex(nxp+15),prey(nyp*2+15),prextm(nxp+15)
      real wbr(nxp/2+1,nyp),wbi(nxp/2+1,nyp)
      real pxz(nxp+2,nzd),w(nxp+2,nzd),bucorr(nxp,nyp,3)
      real bstart,blength,fstart,fend,dstar,rlam,spanv,re
      real u0low,w0low,x0,xs,xlold,xl
      integer fltype,nxin,nxtmp
      real ybl,etabl,df1
      real eta(nyp),fbla(mbla,7+3*scalar)
      real dybla
      integer nbla
      real m1,vsuc
      logical suction,asbl

      integer xbsta,xfend,nxd,x,y,z,i,nzpad

c
c     MPI
c
      integer my_node
c
c     Functions
c
      real, external :: cubip


      call vrffti(nxtmp,prextm,0)
      xbsta=int(real(nxp)*(fstart+.5*xl)/xl)+1
      xfend=int(real(nxp)*(fend  +.5*xl)/xl)+1
      nxin=nxtmp*2/(2+nfxd)
      nxd=(nxp-nxtmp)/2
      if (my_node.eq.0) then
         write(*,*)'Spatial Box Extension is performed.'
         write(*,*) 'New box length is',xl/dstar
      end if
      if (abs(1.-xl/xlold*real(nxin)/real(nx)).gt.1.e-3
     &     .or.xl.lt.xlold) then
         if (my_node.eq.0) then
            write(*,*)'The new box size must be bigger than '
     &           //'the original.'
            write(*,*)'Or new dx has to be same as the old one.'
         end if
         call stopnow(676096)
      end if
      do i=1,3
         do y=1,nyp
            do x=1,nxp
               bucorr(x,y,i)=0.0
            end do
         end do
      end do

      if (fltype.eq.9) then
         if (xlold+fstart.lt.xl*.5) then
            ybl=1.+eta(y)
            etabl=ybl*sqrt(re/(2.*x0))
            df1=cubip(etabl,fbla(1,2),dybla,nbla)
            do y=1,nyp
               do x=nxtmp+xbsta,nxp
                  bucorr(x,y,1)=df1+u0low
               end do
               do x=1,xbsta-1
                  bucorr(x,y,1)=df1+u0low
               end do
            end do
         else
            do y=1,nyp
               ybl=1.+eta(y)
               etabl=ybl*sqrt(re/(2.*x0))
               df1=cubip(etabl,fbla(1,2),dybla,nbla)
               do x=xbsta-2*nxd,xbsta-1
                  bucorr(x,y,1)=df1+u0low
               end do
            end do
         end if
      else
c
c     If fltype not 9 we need a baseflow for both the new and
c     the old field size
c
         call cbflow(bu1,bu2,xs,x0,eta,
     &        fbla,rlam,m1,spanv,re,!!!!nxtmp!!!!,xlold,dstar,
     &        bstart,blength,tabfre,namfre,u0low,w0low,prextm,prey,
     &        wbr,wbi,fltype,my_node,suction,asbl,vsuc,nbla,dybla)

         if (mod(xbsta-nxd,2).eq.1) then
            do i=1,3
               do y=1,nyp
                  do x=xbsta,xfend-1,2
                     bucorr(x,y,i)=-bu1((x-nxd+1)/2,y,i)
                     bucorr(x+1,y,i)=-bu2((x-nxd+1)/2,y,i)
                  end do
                  if (mod(xbsta-xfend,2).eq.0)
     &                 bucorr(xfend,y,i)=-bu1((xfend-nxd+1)/2,y,i)
               end do
            end do
         else
            do i=1,3
               do y=1,nyp
                  do x=xbsta,xfend-1,2
                     bucorr(x,y,i)=-bu2((x-nxd)/2,y,i)
                     bucorr(x+1,y,i)=-bu1((x-nxd)/2+1,y,i)
                  end do
                  if (mod(xbsta-xfend,2).eq.0)
     &                 bucorr(xfend,y,i)=-bu2((xfend-nxd)/2,y,i)
               end do
            end do
         end if

         call cbflow(bu1,bu2,xs,x0,eta,
     &        fbla,rlam,m1,spanv,re,!!!nxp!!!,xl,dstar,
     &        bstart,blength,tabfre,namfre,u0low,w0low,prex,prey,
     &        wbr,wbi,fltype,my_node,suction,asbl,vsuc,nbla,dybla)

         do i=1,3
            do y=1,nyp
               do x=1,nxp/2
                  bucorr(2*x-1,y,i)=bucorr(2*x-1,y,i)+bu1(x,y,i)
                  bucorr(2*x,y,i)=bucorr(2*x,y,i)+bu2(x,y,i)
               end do
            end do
         end do
      end if

      nzpad=(nzp-nz)

      do i=1,3
         sym=i.le.2
         do y=1,nyp
c
c     Get plane form core
c
            do z=1,nz/2
               do x=1,nxin/2
                  pxz(2*x-1,z)=ur(x,y,z,i)
                  pxz(2*x,z)=ui(x,y,z,i)
               end do
            end do
            if (nfzsym.eq.0) then
               do z=nz/2+1,nz
                  do x=1,nxin/2
                     pxz(2*x-1,z+nzpad)=ur(x,y,z,i)
                     pxz(2*x,z+nzpad)=ui(x,y,z,i)
                  end do
               end do
            end if
c
c     Pad with zeros
c
            do z=(nz+1)/2+1,min(nzpc,nzp+1-nz/2)
               do x=1,nxp+2
                  pxz(x,z)=0.0
               end do
            end do
            do z=1,nzpc
               do x=nxin+1,nxp+2
                  pxz(x,z)=0.0
               end do
            end do
c
c     Transform to physical space
c     first complex transform in z direction
c
            if (nfzsym.eq.0) then
               call vcfftb(pxz,pxz(2,1),w,w(2,1),nzp,nxtmp/2,
     &              nxp+2,2,prez)
            else
               if (sym) then
                  call vcffts(pxz,pxz(2,1),w,nzst,nxtmp/2,nxp+2,2,pres)
               else
                  call vcftab(pxz(1,2),pxz(2,2),w,w(2,1),nzat,nxtmp/2,
     &                 nxp+2,2,prea)
               end if
            end if
c
c     Then half complex to real transform in x-direction
c
            call vrfftb(pxz,pxz(2,1),w,w(2,1),nxtmp,nzpc,2,nxp+2,prextm)
c
c     Extend box in streamwise direction
c
            if (xlold+fstart.le.xl*.5) then
               do z=1,nzpc
                  do x=nxp,xbsta+nxtmp,-1
                     pxz(x,z)=bucorr(x,y,i)
                  end do
                  do x=xbsta+nxtmp-1,nxtmp+nxd+1,-1
                     pxz(x,z)=pxz(x-nxtmp-nxd,z)
                  end do
                  do x=nxtmp+nxd,xfend+1,-1
                     pxz(x,z)=pxz(x-nxd,z)
                  end do
                  do x=xfend,xbsta,-1
                     pxz(x,z)=pxz(x-nxd,z)+bucorr(x,y,i)
                  end do
                  do x=xbsta-1,1,-1
                     pxz(x,z)=bucorr(x,y,i)
                  end do
               end do
            else
               do z=1,nzpc
                  do x=nxp,nxtmp+nxd+1,-1
                     pxz(x,z)=pxz(x-nxtmp-nxd,z)
                  end do
                  do x=nxtmp+nxd,xfend+1,-1
                     pxz(x,z)=pxz(x-nxd,z)
                  end do
                  do x=xfend,xbsta,-1
                     pxz(x,z)=pxz(x-nxd,z)+bucorr(x,y,i)
                  end do
                  do x=1,xbsta-2*nxd-1
                     pxz(x,z)=pxz(x+nxd,z)
                  end do
                  do x=xbsta-2*nxd,xbsta-1
                     pxz(x,z)=bucorr(x,y,i)
                  end do
               end do
            end if
c
c     Transform back to spectral space
c     real to half complex transform in x-direction first
c
            call vrfftf(pxz,pxz(2,1),w,w(2,1),nxp,nzpc,2,nxp+2,prex)
c
c     Then complex transform in z direction
c
            if (nfzsym.eq.0) then
               call vcfftf(pxz,pxz(2,1),w,w(2,1),nzp,nxp/2,nxp+2,2,prez)
            else
               if (sym) then
                  call vcffts(pxz,pxz(2,1),w,nzst,nxp/2,nxp+2,2,pres)
               else
                  call vcftaf(pxz(1,2),pxz(2,2),w,w(2,1),nzat,nxp/2,
     &                 nxp+2,2,prea)
               end if
            end if
c
c     Normalize
c
            do z=1,nzd
               do x=1,nxp+2
                  pxz(x,z)=pxz(x,z)/(real(nxp)*real(nzp))
               end do
            end do
c
c     Put plane in core
c
            do z=1,nz/2
               do x=1,nx/2
                  ur(x,y,z,i)=pxz(2*x-1,z)
                  ui(x,y,z,i)=pxz(2*x,z)
               end do
            end do
            if (nfzsym.eq.0) then
c
c     Remove the new oddball modes
c
               if (mod(nzc,2).eq.0) then
                  do x=1,nx/2
                     pxz(2*x-1,nz/2+1+nzpad)=0.
                     pxz(2*x,nz/2+1+nzpad)=0.
                  end do
               end if
               do z=nz/2+1,nz
                  do x=1,nx/2
                     ur(x,y,z,i)=pxz(2*x-1,z+nzpad)
                     ui(x,y,z,i)=pxz(2*x,z+nzpad)
                  end do
               end do
            end if
         end do
      end do

      end subroutine lenbox
