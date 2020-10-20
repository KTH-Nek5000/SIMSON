c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
c     pi                               = 3.1415926535897932385
c     taur                             =-0.5
c     taui = sqrt(0.75)                = 0.86602540378443864676
c     tr11 = sin(4.*atan(1.)/10.)      = 0.30901699437494742410
c     ti11 = sin(2.* 4.*atan(1.)/5.)   = 0.95105651629515357212
c     tr12 = -sin(4.*atan(1.)/10.)-0.5 =-0.80901699437494742410
c     ti12 = sin(4.*atan(1.)/5.)       = 0.58778525229247312917
c     tr13 = sin(4.*atan(1.)/10.)+0.25 = 0.55901699437494742410
c
c     For info see vecfft.txt
c
      subroutine fft_identify(string)
      implicit none
      character(len=*) string

      string = 'cvecfft_orig'

      end subroutine fft_identify

      subroutine vcfftb (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n
      real cr(1+inc*(n-1)+jmp*(m-1)),ci(1+inc*(n-1)+jmp*(m-1))
      real wr(1+inc*(n-1)+jmp*(m-1)),wi(1+inc*(n-1)+jmp*(m-1))
      real wsave(15+2*n)
      if (n .eq. 1) return
      call cfftb1 (cr,ci,wr,wi,m,n,inc,jmp,wsave,wsave(2*n+1))
      return
      end

      subroutine cfftb1 (cr,ci,wr,wi,m,n,inc,jmp,wa,ifac)
      implicit none
      integer inc,jmp,m,n,ifac(15)
      real cr(1+inc*(n-1)+jmp*(m-1)),ci(1+inc*(n-1)+jmp*(m-1))
      real wr(1+inc*(n-1)+jmp*(m-1)),wi(1+inc*(n-1)+jmp*(m-1))
      real wa(2*n)
      integer nf,na,l1,k1,ido,l2,idot,idl1,ip
      integer iw,ix2,ix3,ix4,ix5,ix6,ix7
      integer i,j
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         ix5 = ix4+idot
         ix6 = ix5+idot
         ix7 = ix6+idot
         if (ip.gt.9.or.ip.eq.7) then
           write(6,*) 'factor',ip,'not implemented'
           stop
         end if
         if (na .eq. 0) then
           if (ip .eq. 4) call passb4 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3))
           if (ip .eq. 6) call passb6 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5))
           if (ip .eq. 5) call passb5 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4))
           if (ip .eq. 8) call passb8 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &         wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5),wa(ix6),wa(ix7))
           if (ip .eq. 3) call passb3 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2))
           if (ip .eq. 2) call passb2 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw))
         else
           if (ip .eq. 4) call passb4 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3))
           if (ip .eq. 6) call passb6 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5))
           if (ip .eq. 5) call passb5 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4))
           if (ip .eq. 8) call passb8 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &         wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5),wa(ix6),wa(ix7))
           if (ip .eq. 3) call passb3 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2))
           if (ip .eq. 2) call passb2 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw))
         end if
         na=1-na
         l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      do 117 i=0,n-1
      do 117 j=0,m-1
         cr(1+inc*i+jmp*j) = wr(1+inc*i+jmp*j)
         ci(1+inc*i+jmp*j) = wi(1+inc*i+jmp*j)
  117 continue
      return
      end

      subroutine passb2 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido)
      integer i,j,k,ii,iu
      real tr2,ti2
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*2*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          chr(iu,1) = ccr(ii,1)+ccr(ii,2)
          chr(iu,2) = ccr(ii,1)-ccr(ii,2)
          chi(iu,1) = cci(ii,1)+cci(ii,2)
          chi(iu,2) = cci(ii,1)-cci(ii,2)
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*2*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          chr(iu,1) = ccr(ii,1)+ccr(ii,2)
          tr2 = ccr(ii,1)-ccr(ii,2)
          chi(iu,1) = cci(ii,1)+cci(ii,2)
          ti2 = cci(ii,1)-cci(ii,2)
          chi(iu,2) = wa1(2*i-1)*ti2+wa1(2*i)*tr2
          chr(iu,2) = wa1(2*i-1)*tr2-wa1(2*i)*ti2
  103   continue
      end if
      return
      end

      subroutine passb3 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1,wa2)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido)
      integer i,j,k,ii,iu
      real taur,taui
      real tr2,ti2,cr2,ci2,cr3,ci3,dr2,di2,dr3,di3
      parameter (taur=-.5,taui=.866025403784439)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*3*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          tr2 = ccr(ii,2)+ccr(ii,3)
          cr2 = ccr(ii,1)+taur*tr2
          chr(iu,1) = ccr(ii,1)+tr2
          ti2 = cci(ii,2)+cci(ii,3)
          ci2 = cci(ii,1)+taur*ti2
          chi(iu,1) = cci(ii,1)+ti2
          cr3 = taui*(ccr(ii,2)-ccr(ii,3))
          ci3 = taui*(cci(ii,2)-cci(ii,3))
          chr(iu,2) = cr2-ci3
          chr(iu,3) = cr2+ci3
          chi(iu,2) = ci2+cr3
          chi(iu,3) = ci2-cr3
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*3*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr2 = ccr(ii,2)+ccr(ii,3)
          cr2 = ccr(ii,1)+taur*tr2
          chr(iu,1) = ccr(ii,1)+tr2
          ti2 = cci(ii,2)+cci(ii,3)
          ci2 = cci(ii,1)+taur*ti2
          chi(iu,1) = cci(ii,1)+ti2
          cr3 = taui*(ccr(ii,2)-ccr(ii,3))
          ci3 = taui*(cci(ii,2)-cci(ii,3))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          chi(iu,2) = wa1(2*i-1)*di2+wa1(2*i)*dr2
          chr(iu,2) = wa1(2*i-1)*dr2-wa1(2*i)*di2
          chi(iu,3) = wa2(2*i-1)*di3+wa2(2*i)*dr3
          chr(iu,3) = wa2(2*i-1)*dr3-wa2(2*i)*di3
  103   continue
      end if
      return
      end

      subroutine passb4 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1,wa2,wa3)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido)
      integer i,j,k,ii,iu
      real tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,cr3,ci3
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*4*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          ti1 = cci(ii,1)-cci(ii,3)
          ti2 = cci(ii,1)+cci(ii,3)
          ti3 = cci(ii,2)+cci(ii,4)
          tr4 = cci(ii,4)-cci(ii,2)
          chi(iu,1) = ti2+ti3
          chi(iu,3) = ti2-ti3
          chr(iu,2) = ccr(ii,1)-ccr(ii,3)+tr4
          chr(iu,4) = ccr(ii,1)-ccr(ii,3)-tr4
C          tr1 = ccr(ii,1)-ccr(ii,3)
C          tr2 = ccr(ii,1)+ccr(ii,3)
C          ti4 = ccr(ii,2)-ccr(ii,4)
C          tr3 = ccr(ii,2)+ccr(ii,4)
          chi(iu,2) = ti1+(ccr(ii,2)-ccr(ii,4))
          chi(iu,4) = ti1-(ccr(ii,2)-ccr(ii,4))
          chr(iu,1) = ccr(ii,1)+ccr(ii,3)+(ccr(ii,2)+ccr(ii,4))
          chr(iu,3) = ccr(ii,1)+ccr(ii,3)-(ccr(ii,2)+ccr(ii,4))
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*4*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr1 = ccr(ii,1)-ccr(ii,3)
          tr2 = ccr(ii,1)+ccr(ii,3)
          tr3 = ccr(ii,2)+ccr(ii,4)
          cr3 = tr2-tr3
          ti4 = ccr(ii,2)-ccr(ii,4)
          chr(iu,1) = tr2+tr3
          ti1 = cci(ii,1)-cci(ii,3)
          ti2 = cci(ii,1)+cci(ii,3)
          ti3 = cci(ii,2)+cci(ii,4)
          tr4 = cci(ii,4)-cci(ii,2)
          chi(iu,1) = ti2+ti3
          ci3 = ti2-ti3
          chr(iu,3) = wa2(2*i-1)*cr3-wa2(2*i)*ci3
          chi(iu,3) = wa2(2*i-1)*ci3+wa2(2*i)*cr3
C          cr2 = tr1+tr4
C          cr4 = tr1-tr4
C          ci2 = ti1+ti4
C          ci4 = ti1-ti4
          chr(iu,2) = wa1(2*i-1)*(tr1+tr4)-wa1(2*i)*(ti1+ti4)
          chi(iu,2) = wa1(2*i-1)*(ti1+ti4)+wa1(2*i)*(tr1+tr4)
          chr(iu,4) = wa3(2*i-1)*(tr1-tr4)-wa3(2*i)*(ti1-ti4)
          chi(iu,4) = wa3(2*i-1)*(ti1-ti4)+wa3(2*i)*(tr1-tr4)
  103   continue
      end if
      return
      end

      subroutine passb5(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido)
      integer i,j,k,ii,iu
      real tr11,ti11,tr12,ti12
      real tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5,cr2,ci2,cr3,ci3,cr4,ci4
      real cr5,ci5

      parameter(tr11= .309016994374947)
      parameter(ti11= .951056516295154)
      parameter(tr12=-.809016994374947)
      parameter(ti12= .587785252292473)

      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*5*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          ti5 = cci(ii,2)-cci(ii,5)
          ti2 = cci(ii,2)+cci(ii,5)
          ti4 = cci(ii,3)-cci(ii,4)
          ti3 = cci(ii,3)+cci(ii,4)
          chi(iu,1) = cci(ii,1)+ti2+ti3
          tr5 = ccr(ii,2)-ccr(ii,5)
          tr2 = ccr(ii,2)+ccr(ii,5)
          tr4 = ccr(ii,3)-ccr(ii,4)
          tr3 = ccr(ii,3)+ccr(ii,4)
          chr(iu,1) = ccr(ii,1)+tr2+tr3
          cr2 = ccr(ii,1)+tr11*tr2+tr12*tr3
          cr3 = ccr(ii,1)+tr12*tr2+tr11*tr3
          ci2 = cci(ii,1)+tr11*ti2+tr12*ti3
          ci3 = cci(ii,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          cr4 = ti12*tr5-ti11*tr4
          ci5 = ti11*ti5+ti12*ti4
          ci4 = ti12*ti5-ti11*ti4
          chr(iu,2) = cr2-ci5
          chr(iu,5) = cr2+ci5
          chi(iu,2) = ci2+cr5
          chi(iu,5) = ci2-cr5
          chr(iu,3) = cr3-ci4
          chr(iu,4) = cr3+ci4
          chi(iu,3) = ci3+cr4
          chi(iu,4) = ci3-cr4
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*5*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          ti5 = cci(ii,2)-cci(ii,5)
          ti2 = cci(ii,2)+cci(ii,5)
          ti4 = cci(ii,3)-cci(ii,4)
          ti3 = cci(ii,3)+cci(ii,4)
          tr5 = ccr(ii,2)-ccr(ii,5)
          tr2 = ccr(ii,2)+ccr(ii,5)
          tr4 = ccr(ii,3)-ccr(ii,4)
          tr3 = ccr(ii,3)+ccr(ii,4)
          chr(iu,1) = ccr(ii,1)+tr2+tr3
          chi(iu,1) = cci(ii,1)+ti2+ti3
          cr2 = ccr(ii,1)+tr11*tr2+tr12*tr3
          cr3 = ccr(ii,1)+tr12*tr2+tr11*tr3
          ci2 = cci(ii,1)+tr11*ti2+tr12*ti3
          ci3 = cci(ii,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          cr4 = ti12*tr5-ti11*tr4
          ci5 = ti11*ti5+ti12*ti4
          ci4 = ti12*ti5-ti11*ti4
C          dr3 = cr3-ci4
C          dr4 = cr3+ci4
C          di3 = ci3+cr4
C          di4 = ci3-cr4
          chr(iu,3) = wa2(2*i-1)*(cr3-ci4)-wa2(2*i)*(ci3+cr4)
          chi(iu,3) = wa2(2*i-1)*(ci3+cr4)+wa2(2*i)*(cr3-ci4)
          chr(iu,4) = wa3(2*i-1)*(cr3+ci4)-wa3(2*i)*(ci3-cr4)
          chi(iu,4) = wa3(2*i-1)*(ci3-cr4)+wa3(2*i)*(cr3+ci4)
C          dr5 = cr2+ci5
C          dr2 = cr2-ci5
C          di5 = ci2-cr5
C          di2 = ci2+cr5
          chr(iu,2) = wa1(2*i-1)*(cr2-ci5)-wa1(2*i)*(ci2+cr5)
          chi(iu,2) = wa1(2*i-1)*(ci2+cr5)+wa1(2*i)*(cr2-ci5)
          chr(iu,5) = wa4(2*i-1)*(cr2+ci5)-wa4(2*i)*(ci2-cr5)
          chi(iu,5) = wa4(2*i-1)*(ci2-cr5)+wa4(2*i)*(cr2+ci5)
  103   continue
      end if
      return
      end

      subroutine passb6(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4,wa5)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido),wa5(2*ido)
      integer i,j,k,ii,iu
      real taur,taui
      real tr2,ti2,tr5,ti5
      real cr2,ci2,cr3,ci3,cr5,ci5,cr6,ci6
      real dr1,di1,dr2,di2,dr3,di3,dr4,di4,dr5,di5,dr6,di6
      real er2,ei2,er3,ei3,er4,ei4,er5,ei5,er6,ei6
      parameter (taur=-.5,taui=.866025403784439)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*6*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          tr2 = ccr(ii,3)+ccr(ii,5)
          cr3 = taui*(ccr(ii,3)-ccr(ii,5))
          cr2 = ccr(ii,1)+taur*tr2
          dr1 = ccr(ii,1)+tr2
          ti2 = cci(ii,3)+cci(ii,5)
          ci3 = taui*(cci(ii,3)-cci(ii,5))
          ci2 = cci(ii,1)+taur*ti2
          di1 = cci(ii,1)+ti2
C          dr3 = cr2-ci3
C          dr5 = cr2+ci3
          tr5 = ccr(ii,6)+ccr(ii,2)
          cr6 = taui*(ccr(ii,6)-ccr(ii,2))
          cr5 = ccr(ii,4)+taur*tr5
          dr2 = ccr(ii,4)+tr5
          ti5 = cci(ii,6)+cci(ii,2)
          ci6 = taui*(cci(ii,6)-cci(ii,2))
          ci5 = cci(ii,4)+taur*ti5
          di2 = cci(ii,4)+ti5
          chi(iu,1) = di1+di2
          chi(iu,4) = di1-di2
          chr(iu,1) = dr1+dr2
          chr(iu,4) = dr1-dr2
C          dr4 = cr5-ci6
C          dr6 = cr5+ci6
C          di3 = ci2+cr3
C          di5 = ci2-cr3
C          di4 = ci5+cr6
C          di6 = ci5-cr6
          chi(iu,5) = ci2+cr3+(ci5+cr6)
          chi(iu,2) = ci2+cr3-(ci5+cr6)
          chi(iu,3) = ci2-cr3+(ci5-cr6)
          chi(iu,6) = ci2-cr3-(ci5-cr6)
          chr(iu,5) = cr2-ci3+(cr5-ci6)
          chr(iu,2) = cr2-ci3-(cr5-ci6)
          chr(iu,3) = cr2+ci3+(cr5+ci6)
          chr(iu,6) = cr2+ci3-(cr5+ci6)
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*6*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr2 = ccr(ii,3)+ccr(ii,5)
          cr3 = taui*(ccr(ii,3)-ccr(ii,5))
          cr2 = ccr(ii,1)+taur*tr2
          dr1 = ccr(ii,1)+tr2
          ti2 = cci(ii,3)+cci(ii,5)
          ci3 = taui*(cci(ii,3)-cci(ii,5))
          ci2 = cci(ii,1)+taur*ti2
          di1 = cci(ii,1)+ti2
          dr3 = cr2-ci3
          dr5 = cr2+ci3
          di3 = ci2+cr3
          di5 = ci2-cr3
          tr5 = ccr(ii,6)+ccr(ii,2)
          cr6 = taui*(ccr(ii,6)-ccr(ii,2))
          cr5 = ccr(ii,4)+taur*tr5
          dr2 = ccr(ii,4)+tr5
          ti5 = cci(ii,6)+cci(ii,2)
          ci6 = taui*(cci(ii,6)-cci(ii,2))
          ci5 = cci(ii,4)+taur*ti5
          di2 = cci(ii,4)+ti5
          chr(iu,1) = dr1+dr2
          er4 = dr1-dr2
          chi(iu,1) = di1+di2
          ei4 = di1-di2
          dr4 = cr5-ci6
          dr6 = cr5+ci6
          di4 = ci5+cr6
          di6 = ci5-cr6
          chr(iu,4) = wa3(2*i-1)*er4-wa3(2*i)*ei4
          chi(iu,4) = wa3(2*i-1)*ei4+wa3(2*i)*er4
          er5 = dr3+dr4
          er2 = dr3-dr4
          ei5 = di3+di4
          ei2 = di3-di4
          chr(iu,2) = wa1(2*i-1)*er2-wa1(2*i)*ei2
          chi(iu,2) = wa1(2*i-1)*ei2+wa1(2*i)*er2
          chr(iu,5) = wa4(2*i-1)*er5-wa4(2*i)*ei5
          chi(iu,5) = wa4(2*i-1)*ei5+wa4(2*i)*er5
          er3 = dr5+dr6
          er6 = dr5-dr6
          ei3 = di5+di6
          ei6 = di5-di6
          chr(iu,3) = wa2(2*i-1)*er3-wa2(2*i)*ei3
          chi(iu,3) = wa2(2*i-1)*ei3+wa2(2*i)*er3
          chr(iu,6) = wa5(2*i-1)*er6-wa5(2*i)*ei6
          chi(iu,6) = wa5(2*i-1)*ei6+wa5(2*i)*er6
  103   continue
      end if
      return
      end

      subroutine passb8(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4,wa5,wa6,wa7)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido),wa5(2*ido)
      real wa6(2*ido),wa7(2*ido)
      integer i,j,k,ii,iu
      real w
      real tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4
      real tr5,ti5,tr6,ti6,tr7,ti7,tr8,ti8
      real sr1,si1,sr2,si2,sr3,si3,sr4,si4
      real sr5,si5,sr6,si6,sr7,si7,sr8,si8
      real qr1,qi1,qr2,qi2,cr2,ci2,cr3,ci3,cr4,ci4
      real cr5,ci5,cr6,ci6,cr7,ci7,cr8,ci8
      parameter (w=.7071067811865476)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*8*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          sr1 = ccr(ii,1)+ccr(ii,5)
          sr2 = ccr(ii,1)-ccr(ii,5)
          sr3 = ccr(ii,2)+ccr(ii,6)
          sr4 = ccr(ii,2)-ccr(ii,6)
          sr5 = ccr(ii,3)+ccr(ii,7)
          sr6 = ccr(ii,3)-ccr(ii,7)
          sr7 = ccr(ii,4)+ccr(ii,8)
          sr8 = ccr(ii,4)-ccr(ii,8)
          tr1 = sr1+sr5
          tr2 = sr1-sr5
          tr3 = sr3+sr7
          tr4 = sr3-sr7
          chr(iu,1) = tr1+tr3
          chr(iu,5) = tr1-tr3
          si1 = cci(ii,1)+cci(ii,5)
          si2 = cci(ii,1)-cci(ii,5)
          si3 = cci(ii,2)+cci(ii,6)
          si4 = cci(ii,2)-cci(ii,6)
          si5 = cci(ii,3)+cci(ii,7)
          si6 = cci(ii,3)-cci(ii,7)
          si7 = cci(ii,4)+cci(ii,8)
          si8 = cci(ii,4)-cci(ii,8)
          ti1 = si1+si5
          ti2 = si1-si5
          ti3 = si3+si7
          ti4 = si3-si7
          chi(iu,1) = ti1+ti3
          chi(iu,5) = ti1-ti3
          chr(iu,3) = tr2-ti4
          chr(iu,7) = tr2+ti4
          chi(iu,3) = ti2+tr4
          chi(iu,7) = ti2-tr4
          tr5 = sr2-si6
          tr6 = sr2+si6
          tr7 = w*(sr4-si4)
          ti7 = w*(si4+sr4)
          tr8 = w*(sr8-si8)
          ti8 = w*(si8+sr8)
          ti5 = si2+sr6
          ti6 = si2-sr6
          qi1 = ti7+tr8
          qr2 = tr8-ti7
          qr1 = tr7-ti8
          qi2 = ti8+tr7
          chr(iu,2) = tr5+qr1
          chr(iu,6) = tr5-qr1
          chi(iu,2) = ti5+qi1
          chi(iu,6) = ti5-qi1
          chr(iu,4) = tr6+qr2
          chr(iu,8) = tr6-qr2
          chi(iu,4) = ti6+qi2
          chi(iu,8) = ti6-qi2
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*8*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          sr1 = ccr(ii,1)+ccr(ii,5)
          sr2 = ccr(ii,1)-ccr(ii,5)
          sr3 = ccr(ii,2)+ccr(ii,6)
          sr4 = ccr(ii,2)-ccr(ii,6)
          sr5 = ccr(ii,3)+ccr(ii,7)
          sr6 = ccr(ii,3)-ccr(ii,7)
          sr7 = ccr(ii,4)+ccr(ii,8)
          sr8 = ccr(ii,4)-ccr(ii,8)
          tr1 = sr1+sr5
          tr2 = sr1-sr5
          tr3 = sr3+sr7
          tr4 = sr3-sr7
          chr(iu,1) = tr1+tr3
          cr5 = tr1-tr3
          si1 = cci(ii,1)+cci(ii,5)
          si2 = cci(ii,1)-cci(ii,5)
          si3 = cci(ii,2)+cci(ii,6)
          si4 = cci(ii,2)-cci(ii,6)
          si5 = cci(ii,3)+cci(ii,7)
          si6 = cci(ii,3)-cci(ii,7)
          si7 = cci(ii,4)+cci(ii,8)
          si8 = cci(ii,4)-cci(ii,8)
          ti1 = si1+si5
          ti2 = si1-si5
          ti3 = si3+si7
          ti4 = si3-si7
          chi(iu,1) = ti1+ti3
          ci5 = ti1-ti3
          chr(iu,5) = wa4(2*i-1)*cr5-wa4(2*i)*ci5
          chi(iu,5) = wa4(2*i-1)*ci5+wa4(2*i)*cr5
          cr3 = tr2-ti4
          cr7 = tr2+ti4
          ci3 = ti2+tr4
          ci7 = ti2-tr4
          chr(iu,3) = wa2(2*i-1)*cr3-wa2(2*i)*ci3
          chi(iu,3) = wa2(2*i-1)*ci3+wa2(2*i)*cr3
          chr(iu,7) = wa6(2*i-1)*cr7-wa6(2*i)*ci7
          chi(iu,7) = wa6(2*i-1)*ci7+wa6(2*i)*cr7
          tr5 = sr2-si6
          tr6 = sr2+si6
          tr7 = w*(sr4-si4)
          ti7 = w*(si4+sr4)
          tr8 = w*(sr8-si8)
          ti8 = w*(si8+sr8)
          ti5 = si2+sr6
          ti6 = si2-sr6
          qi1 = ti7+tr8
          qr2 = tr8-ti7
          qr1 = tr7-ti8
          qi2 = ti8+tr7
          cr2 = tr5+qr1
          cr6 = tr5-qr1
          ci2 = ti5+qi1
          ci6 = ti5-qi1
          chr(iu,2) = wa1(2*i-1)*cr2-wa1(2*i)*ci2
          chi(iu,2) = wa1(2*i-1)*ci2+wa1(2*i)*cr2
          chr(iu,6) = wa5(2*i-1)*cr6-wa5(2*i)*ci6
          chi(iu,6) = wa5(2*i-1)*ci6+wa5(2*i)*cr6
          cr4 = tr6+qr2
          cr8 = tr6-qr2
          ci4 = ti6+qi2
          ci8 = ti6-qi2
          chr(iu,4) = wa3(2*i-1)*cr4-wa3(2*i)*ci4
          chi(iu,4) = wa3(2*i-1)*ci4+wa3(2*i)*cr4
          chr(iu,8) = wa7(2*i-1)*cr8-wa7(2*i)*ci8
          chi(iu,8) = wa7(2*i-1)*ci8+wa7(2*i)*cr8
  103   continue
      end if
      return
      end

cvd$g nolstval
      subroutine vcfftf (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n
      real cr(1+inc*(n-1)+jmp*(m-1)),ci(1+inc*(n-1)+jmp*(m-1))
      real wr(1+inc*(n-1)+jmp*(m-1)),wi(1+inc*(n-1)+jmp*(m-1))
      real wsave(15+2*n)
      if (n .eq. 1) return
      call cfftf1 (cr,ci,wr,wi,m,n,inc,jmp,wsave,wsave(2*n+1))
      return
      end

      subroutine cfftf1 (cr,ci,wr,wi,m,n,inc,jmp,wa,ifac)
      implicit none
      integer inc,jmp,m,n,ifac(15),i,j
      real cr(1+inc*(n-1)+jmp*(m-1)),ci(1+inc*(n-1)+jmp*(m-1))
      real wr(1+inc*(n-1)+jmp*(m-1)),wi(1+inc*(n-1)+jmp*(m-1))
      real wa(2*n)
      integer nf,na,l1,k1,ido,l2,idot,idl1,ip
      integer iw,ix2,ix3,ix4,ix5,ix6,ix7
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         ix5 = ix4+idot
         ix6 = ix5+idot
         ix7 = ix6+idot
         if (ip.gt.9.or.ip.eq.7) then
           write(6,*) 'factor',ip,'not implemented'
           stop
         end if
         if (na .eq. 0) then
           if (ip .eq. 4) call passf4 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3))
           if (ip .eq. 6) call passf6 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5))
           if (ip .eq. 5) call passf5 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4))
           if (ip .eq. 8) call passf8 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &         wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5),wa(ix6),wa(ix7))
           if (ip .eq. 3) call passf3 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2))
           if (ip .eq. 2) call passf2 (cr,ci,wr,wi,inc,jmp,m,ido,l1,
     &                    wa(iw))
         else
           if (ip .eq. 4) call passf4 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3))
           if (ip .eq. 6) call passf6 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5))
           if (ip .eq. 5) call passf5 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2),wa(ix3),wa(ix4))
           if (ip .eq. 8) call passf8 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &         wa(iw),wa(ix2),wa(ix3),wa(ix4),wa(ix5),wa(ix6),wa(ix7))
           if (ip .eq. 3) call passf3 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw),wa(ix2))
           if (ip .eq. 2) call passf2 (wr,wi,cr,ci,inc,jmp,m,ido,l1,
     &                    wa(iw))
         end if
         na=1-na
         l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
      do 117 i=0,n-1
      do 117 j=0,m-1
         cr(1+inc*i+jmp*j) = wr(1+inc*i+jmp*j)
         ci(1+inc*i+jmp*j) = wi(1+inc*i+jmp*j)
  117 continue
      return
      end

      subroutine passf2 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido)
      integer i,j,k,ii,iu
      real tr2,ti2
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*2*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          chr(iu,1) = ccr(ii,1)+ccr(ii,2)
          chr(iu,2) = ccr(ii,1)-ccr(ii,2)
          chi(iu,1) = cci(ii,1)+cci(ii,2)
          chi(iu,2) = cci(ii,1)-cci(ii,2)
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*2*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          chr(iu,1) = ccr(ii,1)+ccr(ii,2)
          tr2 = ccr(ii,1)-ccr(ii,2)
          chi(iu,1) = cci(ii,1)+cci(ii,2)
          ti2 = cci(ii,1)-cci(ii,2)
          chi(iu,2) = wa1(2*i-1)*ti2-wa1(2*i)*tr2
          chr(iu,2) = wa1(2*i-1)*tr2+wa1(2*i)*ti2
  103   continue
      end if
      return
      end

      subroutine passf3 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1,wa2)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido)
      integer i,j,k,ii,iu
      real taur,taui
      real tr2,ti2,cr2,ci2,cr3,ci3,dr2,di2,dr3,di3
      parameter (taur=-.5,taui=-.866025403784439)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*3*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          tr2 = ccr(ii,2)+ccr(ii,3)
          cr2 = ccr(ii,1)+taur*tr2
          chr(iu,1) = ccr(ii,1)+tr2
          ti2 = cci(ii,2)+cci(ii,3)
          ci2 = cci(ii,1)+taur*ti2
          chi(iu,1) = cci(ii,1)+ti2
          cr3 = taui*(ccr(ii,2)-ccr(ii,3))
          ci3 = taui*(cci(ii,2)-cci(ii,3))
          chr(iu,2) = cr2-ci3
          chr(iu,3) = cr2+ci3
          chi(iu,2) = ci2+cr3
          chi(iu,3) = ci2-cr3
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*3*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr2 = ccr(ii,2)+ccr(ii,3)
          cr2 = ccr(ii,1)+taur*tr2
          chr(iu,1) = ccr(ii,1)+tr2
          ti2 = cci(ii,2)+cci(ii,3)
          ci2 = cci(ii,1)+taur*ti2
          chi(iu,1) = cci(ii,1)+ti2
          cr3 = taui*(ccr(ii,2)-ccr(ii,3))
          ci3 = taui*(cci(ii,2)-cci(ii,3))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          chi(iu,2) = wa1(2*i-1)*di2-wa1(2*i)*dr2
          chr(iu,2) = wa1(2*i-1)*dr2+wa1(2*i)*di2
          chi(iu,3) = wa2(2*i-1)*di3-wa2(2*i)*dr3
          chr(iu,3) = wa2(2*i-1)*dr3+wa2(2*i)*di3
  103   continue
      end if
      return
      end

      subroutine passf4 (ccr,cci,chr,chi,inc,jmp,m,ido,l1,wa1,wa2,wa3)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido)
      integer i,j,k,ii,iu
      real tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4,cr3,ci3
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*4*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          ti1 = cci(ii,1)-cci(ii,3)
          ti2 = cci(ii,1)+cci(ii,3)
          ti3 = cci(ii,2)+cci(ii,4)
          chi(iu,1) = ti2+ti3
          chi(iu,3) = ti2-ti3
          tr4 = cci(ii,2)-cci(ii,4)
          chr(iu,2) = ccr(ii,1)-ccr(ii,3)+tr4
          chr(iu,4) = ccr(ii,1)-ccr(ii,3)-tr4
C          tr1 = ccr(ii,1)-ccr(ii,3)
C          tr2 = ccr(ii,1)+ccr(ii,3)
C          ti4 = ccr(ii,4)-ccr(ii,2)
C          tr3 = ccr(ii,2)+ccr(ii,4)
          chi(iu,2) = (cci(ii,1)-cci(ii,3))+(ccr(ii,4)-ccr(ii,2))
          chi(iu,4) = (cci(ii,1)-cci(ii,3))-(ccr(ii,4)-ccr(ii,2))
          chr(iu,1) = ccr(ii,1)+ccr(ii,3)+(ccr(ii,2)+ccr(ii,4))
          chr(iu,3) = ccr(ii,1)+ccr(ii,3)-(ccr(ii,2)+ccr(ii,4))
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*4*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr1 = ccr(ii,1)-ccr(ii,3)
          tr2 = ccr(ii,1)+ccr(ii,3)
          tr3 = ccr(ii,2)+ccr(ii,4)
          cr3 = tr2-tr3
          ti4 = ccr(ii,4)-ccr(ii,2)
          chr(iu,1) = tr2+tr3
          ti1 = cci(ii,1)-cci(ii,3)
          ti2 = cci(ii,1)+cci(ii,3)
          ti3 = cci(ii,2)+cci(ii,4)
          tr4 = cci(ii,2)-cci(ii,4)
          chi(iu,1) = ti2+ti3
          ci3 = ti2-ti3
          chr(iu,3) = wa2(2*i-1)*cr3+wa2(2*i)*ci3
          chi(iu,3) = wa2(2*i-1)*ci3-wa2(2*i)*cr3
C          cr2 = tr1+tr4
C          cr4 = tr1-tr4
C          ci2 = ti1+ti4
C          ci4 = ti1-ti4
          chr(iu,2) = wa1(2*i-1)*(tr1+tr4)+wa1(2*i)*(ti1+ti4)
          chi(iu,2) = wa1(2*i-1)*(ti1+ti4)-wa1(2*i)*(tr1+tr4)
          chr(iu,4) = wa3(2*i-1)*(tr1-tr4)+wa3(2*i)*(ti1-ti4)
          chi(iu,4) = wa3(2*i-1)*(ti1-ti4)-wa3(2*i)*(tr1-tr4)
  103   continue
      end if
      return
      end

      subroutine passf5(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido)
      integer i,j,k,ii,iu
      real tr11,ti11,tr12,ti12
      real tr2,ti2,tr3,ti3,tr4,ti4,tr5,ti5,cr2,ci2,cr3,ci3,cr4,ci4
      real cr5,ci5

      parameter(tr11= .309016994374947)
      parameter(ti11=-.951056516295154)
      parameter(tr12=-.809016994374947)
      parameter(ti12=-.587785252292473)

      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*5*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          ti5 = cci(ii,2)-cci(ii,5)
          ti2 = cci(ii,2)+cci(ii,5)
          ti4 = cci(ii,3)-cci(ii,4)
          ti3 = cci(ii,3)+cci(ii,4)
          chi(iu,1) = cci(ii,1)+ti2+ti3
          tr5 = ccr(ii,2)-ccr(ii,5)
          tr2 = ccr(ii,2)+ccr(ii,5)
          tr4 = ccr(ii,3)-ccr(ii,4)
          tr3 = ccr(ii,3)+ccr(ii,4)
          chr(iu,1) = ccr(ii,1)+tr2+tr3
          cr2 = ccr(ii,1)+tr11*tr2+tr12*tr3
          cr3 = ccr(ii,1)+tr12*tr2+tr11*tr3
          ci2 = cci(ii,1)+tr11*ti2+tr12*ti3
          ci3 = cci(ii,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          cr4 = ti12*tr5-ti11*tr4
          ci5 = ti11*ti5+ti12*ti4
          ci4 = ti12*ti5-ti11*ti4
          chr(iu,2) = cr2-ci5
          chr(iu,5) = cr2+ci5
          chi(iu,2) = ci2+cr5
          chi(iu,5) = ci2-cr5
          chr(iu,3) = cr3-ci4
          chr(iu,4) = cr3+ci4
          chi(iu,3) = ci3+cr4
          chi(iu,4) = ci3-cr4
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*5*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          ti5 = cci(ii,2)-cci(ii,5)
          ti2 = cci(ii,2)+cci(ii,5)
          ti4 = cci(ii,3)-cci(ii,4)
          ti3 = cci(ii,3)+cci(ii,4)
          tr5 = ccr(ii,2)-ccr(ii,5)
          tr2 = ccr(ii,2)+ccr(ii,5)
          tr4 = ccr(ii,3)-ccr(ii,4)
          tr3 = ccr(ii,3)+ccr(ii,4)
          chr(iu,1) = ccr(ii,1)+tr2+tr3
          chi(iu,1) = cci(ii,1)+ti2+ti3
          cr2 = ccr(ii,1)+tr11*tr2+tr12*tr3
          cr3 = ccr(ii,1)+tr12*tr2+tr11*tr3
          ci2 = cci(ii,1)+tr11*ti2+tr12*ti3
          ci3 = cci(ii,1)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          cr4 = ti12*tr5-ti11*tr4
          ci5 = ti11*ti5+ti12*ti4
          ci4 = ti12*ti5-ti11*ti4
C          dr3 = cr3-ci4
C          dr4 = cr3+ci4
C          di3 = ci3+cr4
C          di4 = ci3-cr4
          chr(iu,3) = wa2(2*i-1)*(cr3-ci4)+wa2(2*i)*(ci3+cr4)
          chi(iu,3) = wa2(2*i-1)*(ci3+cr4)-wa2(2*i)*(cr3-ci4)
          chr(iu,4) = wa3(2*i-1)*(cr3+ci4)+wa3(2*i)*(ci3-cr4)
          chi(iu,4) = wa3(2*i-1)*(ci3-cr4)-wa3(2*i)*(cr3+ci4)
C          dr5 = cr2+ci5
C          dr2 = cr2-ci5
C          di5 = ci2-cr5
C          di2 = ci2+cr5
          chr(iu,2) = wa1(2*i-1)*(cr2-ci5)+wa1(2*i)*(ci2+cr5)
          chi(iu,2) = wa1(2*i-1)*(ci2+cr5)-wa1(2*i)*(cr2-ci5)
          chr(iu,5) = wa4(2*i-1)*(cr2+ci5)+wa4(2*i)*(ci2-cr5)
          chi(iu,5) = wa4(2*i-1)*(ci2-cr5)-wa4(2*i)*(cr2+ci5)
  103   continue
      end if
      return
      end

      subroutine passf6(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4,wa5)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido),wa5(2*ido)
      integer i,j,k,ii,iu
      real taur,taui
      real tr2,ti2,tr5,ti5
      real cr2,ci2,cr3,ci3,cr5,ci5,cr6,ci6
      real dr1,di1,dr2,di2,dr3,di3,dr4,di4,dr5,di5,dr6,di6
      real er2,ei2,er3,ei3,er4,ei4,er5,ei5,er6,ei6
      parameter (taur=-.5,taui=-.866025403784439)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*6*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          tr2 = ccr(ii,3)+ccr(ii,5)
          cr3 = taui*(ccr(ii,3)-ccr(ii,5))
          cr2 = ccr(ii,1)+taur*tr2
          dr1 = ccr(ii,1)+tr2
          ti2 = cci(ii,3)+cci(ii,5)
          ci3 = taui*(cci(ii,3)-cci(ii,5))
          ci2 = cci(ii,1)+taur*ti2
          di1 = cci(ii,1)+ti2
C          dr3 = cr2-ci3
C          dr5 = cr2+ci3
          tr5 = ccr(ii,6)+ccr(ii,2)
          cr6 = taui*(ccr(ii,6)-ccr(ii,2))
          cr5 = ccr(ii,4)+taur*tr5
          dr2 = ccr(ii,4)+tr5
          ti5 = cci(ii,6)+cci(ii,2)
          ci6 = taui*(cci(ii,6)-cci(ii,2))
          ci5 = cci(ii,4)+taur*ti5
          di2 = cci(ii,4)+ti5
          chi(iu,1) = di1+di2
          chi(iu,4) = di1-di2
          chr(iu,1) = dr1+dr2
          chr(iu,4) = dr1-dr2
C          dr4 = cr5-ci6
C          dr6 = cr5+ci6
C          di3 = ci2+cr3
C          di5 = ci2-cr3
C          di4 = ci5+cr6
C          di6 = ci5-cr6
          chi(iu,5) = ci2+cr3+(ci5+cr6)
          chi(iu,2) = ci2+cr3-(ci5+cr6)
          chi(iu,3) = ci2-cr3+(ci5-cr6)
          chi(iu,6) = ci2-cr3-(ci5-cr6)
          chr(iu,5) = cr2-ci3+(cr5-ci6)
          chr(iu,2) = cr2-ci3-(cr5-ci6)
          chr(iu,3) = cr2+ci3+(cr5+ci6)
          chr(iu,6) = cr2+ci3-(cr5+ci6)
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*6*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          tr2 = ccr(ii,3)+ccr(ii,5)
          cr3 = taui*(ccr(ii,3)-ccr(ii,5))
          cr2 = ccr(ii,1)+taur*tr2
          dr1 = ccr(ii,1)+tr2
          ti2 = cci(ii,3)+cci(ii,5)
          ci3 = taui*(cci(ii,3)-cci(ii,5))
          ci2 = cci(ii,1)+taur*ti2
          di1 = cci(ii,1)+ti2
          dr3 = cr2-ci3
          dr5 = cr2+ci3
          di3 = ci2+cr3
          di5 = ci2-cr3
          tr5 = ccr(ii,6)+ccr(ii,2)
          cr6 = taui*(ccr(ii,6)-ccr(ii,2))
          cr5 = ccr(ii,4)+taur*tr5
          dr2 = ccr(ii,4)+tr5
          ti5 = cci(ii,6)+cci(ii,2)
          ci6 = taui*(cci(ii,6)-cci(ii,2))
          ci5 = cci(ii,4)+taur*ti5
          di2 = cci(ii,4)+ti5
          chr(iu,1) = dr1+dr2
          er4 = dr1-dr2
          chi(iu,1) = di1+di2
          ei4 = di1-di2
          dr4 = cr5-ci6
          dr6 = cr5+ci6
          di4 = ci5+cr6
          di6 = ci5-cr6
          chr(iu,4) = wa3(2*i-1)*er4+wa3(2*i)*ei4
          chi(iu,4) = wa3(2*i-1)*ei4-wa3(2*i)*er4
          er5 = dr3+dr4
          er2 = dr3-dr4
          ei5 = di3+di4
          ei2 = di3-di4
          chr(iu,2) = wa1(2*i-1)*er2+wa1(2*i)*ei2
          chi(iu,2) = wa1(2*i-1)*ei2-wa1(2*i)*er2
          chr(iu,5) = wa4(2*i-1)*er5+wa4(2*i)*ei5
          chi(iu,5) = wa4(2*i-1)*ei5-wa4(2*i)*er5
          er3 = dr5+dr6
          er6 = dr5-dr6
          ei3 = di5+di6
          ei6 = di5-di6
          chr(iu,3) = wa2(2*i-1)*er3+wa2(2*i)*ei3
          chi(iu,3) = wa2(2*i-1)*ei3-wa2(2*i)*er3
          chr(iu,6) = wa5(2*i-1)*er6+wa5(2*i)*ei6
          chi(iu,6) = wa5(2*i-1)*ei6-wa5(2*i)*er6
  103   continue
      end if
      return
      end

      subroutine passf8(ccr,cci,chr,chi,inc,jmp,m,ido,l1,
     & wa1,wa2,wa3,wa4,wa5,wa6,wa7)
      implicit none
      integer inc,jmp,m,ido,l1
      real ccr(inc*ido,1),cci(inc*ido,1)
      real chr(inc*ido*l1,1),chi(inc*ido*l1,1)
      real wa1(2*ido),wa2(2*ido),wa3(2*ido),wa4(2*ido),wa5(2*ido)
      real wa6(2*ido),wa7(2*ido)
      integer i,j,k,ii,iu
      real w
      real tr1,ti1,tr2,ti2,tr3,ti3,tr4,ti4
      real tr5,ti5,tr6,ti6,tr7,ti7,tr8,ti8
      real sr1,si1,sr2,si2,sr3,si3,sr4,si4
      real sr5,si5,sr6,si6,sr7,si7,sr8,si8
      real qr1,qi1,qr2,qi2,cr2,ci2,cr3,ci3,cr4,ci4
      real cr5,ci5,cr6,ci6,cr7,ci7,cr8,ci8
      parameter (w=.7071067811865476)
      if (ido .eq. 1) then
        do 101 k=1,l1
        do 101 j=1,m
          ii=1+jmp*(j-1)+inc*ido*8*(k-1)
          iu=1+jmp*(j-1)+inc*ido*(k-1)
          sr1 = ccr(ii,1)+ccr(ii,5)
          sr2 = ccr(ii,1)-ccr(ii,5)
          sr3 = ccr(ii,2)+ccr(ii,6)
          sr4 = ccr(ii,2)-ccr(ii,6)
          sr5 = ccr(ii,3)+ccr(ii,7)
          sr6 = ccr(ii,3)-ccr(ii,7)
          sr7 = ccr(ii,4)+ccr(ii,8)
          sr8 = ccr(ii,4)-ccr(ii,8)
          tr1 = sr1+sr5
          tr2 = sr1-sr5
          tr3 = sr3+sr7
          tr4 = sr3-sr7
          chr(iu,1) = tr1+tr3
          chr(iu,5) = tr1-tr3
          si1 = cci(ii,1)+cci(ii,5)
          si2 = cci(ii,1)-cci(ii,5)
          si3 = cci(ii,2)+cci(ii,6)
          si4 = cci(ii,2)-cci(ii,6)
          si5 = cci(ii,3)+cci(ii,7)
          si6 = cci(ii,3)-cci(ii,7)
          si7 = cci(ii,4)+cci(ii,8)
          si8 = cci(ii,4)-cci(ii,8)
          ti1 = si1+si5
          ti2 = si1-si5
          ti3 = si3+si7
          ti4 = si3-si7
          chi(iu,1) = ti1+ti3
          chi(iu,5) = ti1-ti3
          chr(iu,3) = tr2+ti4
          chr(iu,7) = tr2-ti4
          chi(iu,3) = ti2-tr4
          chi(iu,7) = ti2+tr4
          tr5 = sr2+si6
          tr6 = sr2-si6
          tr7 = w*(sr4+si4)
          ti7 = w*(si4-sr4)
          tr8 = w*(sr8+si8)
          ti8 = w*(si8-sr8)
          ti5 = si2-sr6
          ti6 = si2+sr6
          qi1 = ti7-tr8
          qr2 = tr8+ti7
          qr1 = tr7+ti8
          qi2 = ti8-tr7
          chr(iu,2) = tr5+qr1
          chr(iu,6) = tr5-qr1
          chi(iu,2) = ti5+qi1
          chi(iu,6) = ti5-qi1
          chr(iu,4) = tr6+qr2
          chr(iu,8) = tr6-qr2
          chi(iu,4) = ti6+qi2
          chi(iu,8) = ti6-qi2
  101   continue
      else
        do 103 k=1,l1
        do 103 i=1,ido
        do 103 j=1,m
          ii=1+jmp*(j-1)+inc*(ido*8*(k-1)+(i-1))
          iu=1+jmp*(j-1)+inc*(ido*(k-1)+(i-1))
          sr1 = ccr(ii,1)+ccr(ii,5)
          sr2 = ccr(ii,1)-ccr(ii,5)
          sr3 = ccr(ii,2)+ccr(ii,6)
          sr4 = ccr(ii,2)-ccr(ii,6)
          sr5 = ccr(ii,3)+ccr(ii,7)
          sr6 = ccr(ii,3)-ccr(ii,7)
          sr7 = ccr(ii,4)+ccr(ii,8)
          sr8 = ccr(ii,4)-ccr(ii,8)
          tr1 = sr1+sr5
          tr2 = sr1-sr5
          tr3 = sr3+sr7
          tr4 = sr3-sr7
          chr(iu,1) = tr1+tr3
          cr5 = tr1-tr3
          si1 = cci(ii,1)+cci(ii,5)
          si2 = cci(ii,1)-cci(ii,5)
          si3 = cci(ii,2)+cci(ii,6)
          si4 = cci(ii,2)-cci(ii,6)
          si5 = cci(ii,3)+cci(ii,7)
          si6 = cci(ii,3)-cci(ii,7)
          si7 = cci(ii,4)+cci(ii,8)
          si8 = cci(ii,4)-cci(ii,8)
          ti1 = si1+si5
          ti2 = si1-si5
          ti3 = si3+si7
          ti4 = si3-si7
          chi(iu,1) = ti1+ti3
          ci5 = ti1-ti3
          chr(iu,5) = wa4(2*i-1)*cr5+wa4(2*i)*ci5
          chi(iu,5) = wa4(2*i-1)*ci5-wa4(2*i)*cr5
          cr3 = tr2+ti4
          cr7 = tr2-ti4
          ci3 = ti2-tr4
          ci7 = ti2+tr4
          chr(iu,3) = wa2(2*i-1)*cr3+wa2(2*i)*ci3
          chi(iu,3) = wa2(2*i-1)*ci3-wa2(2*i)*cr3
          chr(iu,7) = wa6(2*i-1)*cr7+wa6(2*i)*ci7
          chi(iu,7) = wa6(2*i-1)*ci7-wa6(2*i)*cr7
          tr5 = sr2+si6
          tr6 = sr2-si6
          tr7 = w*(sr4+si4)
          ti7 = w*(si4-sr4)
          tr8 = w*(sr8+si8)
          ti8 = w*(si8-sr8)
          ti5 = si2-sr6
          ti6 = si2+sr6
          qi1 = ti7-tr8
          qr2 = tr8+ti7
          qr1 = tr7+ti8
          qi2 = ti8-tr7
          cr2 = tr5+qr1
          cr6 = tr5-qr1
          ci2 = ti5+qi1
          ci6 = ti5-qi1
          chr(iu,2) = wa1(2*i-1)*cr2+wa1(2*i)*ci2
          chi(iu,2) = wa1(2*i-1)*ci2-wa1(2*i)*cr2
          chr(iu,6) = wa5(2*i-1)*cr6+wa5(2*i)*ci6
          chi(iu,6) = wa5(2*i-1)*ci6-wa5(2*i)*cr6
          cr4 = tr6+qr2
          cr8 = tr6-qr2
          ci4 = ti6+qi2
          ci8 = ti6-qi2
          chr(iu,4) = wa3(2*i-1)*cr4+wa3(2*i)*ci4
          chi(iu,4) = wa3(2*i-1)*ci4-wa3(2*i)*cr4
          chr(iu,8) = wa7(2*i-1)*cr8+wa7(2*i)*ci8
          chi(iu,8) = wa7(2*i-1)*ci8-wa7(2*i)*cr8
  103   continue
      end if
      return
      end

      subroutine vcffti (n,wsave,ifail)
      implicit none
      integer ifail,n,ifailn
      integer iw1,iw2
      real wsave(1)
      ifailn=0
      if (n .ne. 1) then
        iw1 = 1
        iw2 = iw1+n+n
        call cffti1 (n,wsave(iw1),wsave(iw2),ifailn)
      end if
      if (ifail.eq.0.and.ifailn.ne.0) then
        write(*,*) 'VCFFTI-ERROR-CANNOT FACTOR',ifailn
        stop
      end if
      if (ifail.ne.0) ifail=ifailn
      return
      end

      subroutine cffti1 (n,wa,ifac,ifail)
      implicit none
      integer n,ifail,nl,nf,j,ntry,np,i,l1,k1,ip,ld,l2,ido
      integer idot,ipm,i1,ii
      real fi
      real wa(1)
      integer ifac(1),ntryh(6),nfac(6)
      data ntryh/2,3,4,5,6,8/
      real tpi,argh,argld,arg
      parameter (tpi = 6.28318530717959)
      ifail=0
      nl = n
      nf = 0
      nfac(6)=0
C first factor n, select large factors first except 8
      do 101 j=5,1,-1
        nfac(j)=0
        ntry=ntryh(j)
 102    if (nl.eq.nl/ntry*ntry) then
          nl=nl/ntry
          nf=nf+1
          nfac(j)=nfac(j)+1
          goto 102
        end if
 101  continue
      if (nl.ne.1) then
        ifail=1
        return
      end if
C make a single factor 8 out of a 4 and a 2 if nf becomes even
C (saves the overhead of copying data at least)
      if (nfac(1).ge.1.and.nfac(3).ge.1.and.nf/2*2.ne.nf) then
        nfac(1)=nfac(1)-1
        nfac(3)=nfac(3)-1
        nfac(6)=nfac(6)+1
        nf=nf-1
      end if
C now insert the factors into ifac but the largest factor last
      np=3
      do 103 j=1,6
        do 104 np=np,np-1+nfac(j)
          ifac(np)=ntryh(j)
 104    continue
 103  continue
      ifac(1) = n
      ifac(2) = nf
C
      argh = tpi/real(n)
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.
            wa(i) = 0.
            ld = ld+l1
            fi=0.
            argld = real(ld)*argh
            do 108 ii=4,idot,2
               i = i+2
               fi=fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
CC            if (ip .le. 5) go to 109
CC            wa(i1-1) = wa(i-1)
CC            wa(i1) = wa(i)
  109    continue
         l1 = l2
  110 continue
      return
      end

cvd$g nolstval
      subroutine vcffts (cr,ci,wr,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n
      real cr(1+inc*(n+1)+jmp*(m-1)),ci(1+inc*(n+1)+jmp*(m-1))
      real wr(1+inc*(n+1)+1+jmp*(m-1))
      real wsave(15+2*n)
      call vcost(cr,wr,n,m,inc,jmp,wsave)
      call vcost(ci,wr,n,m,inc,jmp,wsave)
      return
      end

cvd$g nolstval
      subroutine vcftab (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n,iw1,iw2,np1
      real cr(1+inc*n+jmp*(m-1)),ci(1+inc*n+jmp*(m-1))
      real wr(1+inc*n+jmp*(m-1)),wi(1+inc*n+jmp*(m-1))
      real wsave(15+2*n)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      call cftab1(cr,ci,wr,wi,n,m,inc,jmp,wsave,wsave(iw1),wsave(iw2))
      return
      end

      subroutine cftab1(xr,xi,wr,wi,n,m,inc,jmp,was,war,ifac)
      implicit none
      real war(1),was(1)
      integer ifac(1)
      real xr(inc,1),xi(inc,1),wr(inc,1),wi(inc,1)
      real sqrt3,t1,t2
      integer  n,m,j,jj,jmp,np1,ns2,k,kc,i,inc
      parameter(sqrt3=1.73205080756888)
      if (n.eq.1) then
        do 200 j=1,m
          jj=1+(j-1)*jmp
          t1 = 2.*xr(jj,1)
          xr(jj,1) = -2.*xi(jj,1)
          xi(jj,1) = t1
 200    continue
      else
      np1 = n+1
      ns2 = n/2
      do 205 j=1,m
        jj=1+(j-1)*jmp
        wr(jj,1) = 0.
        wi(jj,1) = 0.
 205  continue
      do 104 k=1,ns2
        kc = np1-k
cdir$ ivdep
        do 202 j=1,m
          jj=1+(j-1)*jmp
          t1 = xr(jj,k)-xr(jj,kc)
          t2 = was(k)*(xr(jj,k)+xr(jj,kc))
          wr(jj,k+1) = t1+t2
          wr(jj,kc+1) = t2-t1
          t1 = xi(jj,k)-xi(jj,kc)
          t2 = was(k)*(xi(jj,k)+xi(jj,kc))
          wi(jj,k+1) = t1+t2
          wi(jj,kc+1) = t2-t1
  202   continue
  104 continue
      do 203 j=1,m
        jj=1+(j-1)*jmp
        wr(jj,ns2+2) = 4.*xr(jj,ns2+1)
        wi(jj,ns2+2) = 4.*xi(jj,ns2+1)
 203  continue
      call rfftf1 (wr,wr(1,2),xr,xr(1,2),np1,m,2*inc,jmp,war,ifac)
      call rfftf1 (wi,wi(1,2),xi,xi(1,2),np1,m,2*inc,jmp,war,ifac)
cdir$ ivdep
      do 204 j=1,m
        jj=1+(j-1)*jmp
        xi(jj,1) = .5*wr(jj,1)
        xr(jj,1) = -.5*wi(jj,1)
        xi(jj,n+1) = 0.
        xr(jj,n+1) = 0.
 204  continue
      do 105 i=3,n,2
cdir$ ivdep
      do 105 j=1,m
         jj=1+(j-1)*jmp
         xi(jj,i-1) = -wr(jj,i+1)
         xi(jj,i) = xi(jj,i-2)+wr(jj,i)
         xr(jj,i-1) = wi(jj,i+1)
         xr(jj,i) = xr(jj,i-2)-wi(jj,i)
  105 continue
      return
      end if
      end

cvd$g nolstval
      subroutine vcftaf (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n,iw1,iw2,np1
      real cr(1+inc*n+jmp*(m-1)),ci(1+inc*n+jmp*(m-1))
      real wr(1+inc*n+jmp*(m-1)),wi(1+inc*n+jmp*(m-1))
      real wsave(15+2*n)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      call cftaf1(cr,ci,wr,wi,n,m,inc,jmp,wsave,wsave(iw1),wsave(iw2))
      return
      end

      subroutine cftaf1(xr,xi,wr,wi,n,m,inc,jmp,was,war,ifac)
      implicit none
      real war(1),was(1)
      integer ifac(1)
      real xr(inc,1),xi(inc,1),wr(inc,1),wi(inc,1)
      real sqrt3
      parameter(sqrt3=1.73205080756888)
      integer n,m,inc,jmp,jj,np1,ns2,k,kc,j,i
      real t1,t2
      if (n.eq.1) then
        do 200 j=1,m
          jj=1+(j-1)*jmp
          t1 = -2.*xr(jj,1)
          xr(jj,1) = 2.*xi(jj,1)
          xi(jj,1) = t1
 200    continue
      else
      np1 = n+1
      ns2 = n/2
      do 205 j=1,m
        jj=1+(j-1)*jmp
        wr(jj,1) = 0.
        wi(jj,1) = 0.
 205  continue
      do 104 k=1,ns2
        kc = np1-k
cdir$ ivdep
        do 202 j=1,m
          jj=1+(j-1)*jmp
          t1 = xr(jj,k)-xr(jj,kc)
          t2 = was(k)*(xr(jj,k)+xr(jj,kc))
          wr(jj,k+1) = t1+t2
          wr(jj,kc+1) = t2-t1
          t1 = xi(jj,k)-xi(jj,kc)
          t2 = was(k)*(xi(jj,k)+xi(jj,kc))
          wi(jj,k+1) = t1+t2
          wi(jj,kc+1) = t2-t1
  202   continue
  104 continue
      do 203 j=1,m
        jj=1+(j-1)*jmp
        wr(jj,ns2+2) = 4.*xr(jj,ns2+1)
        wi(jj,ns2+2) = 4.*xi(jj,ns2+1)
 203  continue
      call rfftf1 (wr,wr(1,2),xr,xr(1,2),np1,m,2*inc,jmp,war,ifac)
      call rfftf1 (wi,wi(1,2),xi,xi(1,2),np1,m,2*inc,jmp,war,ifac)
cdir$ ivdep
      do 204 j=1,m
        jj=1+(j-1)*jmp
        xi(jj,1) = -.5*wr(jj,1)
        xr(jj,1) = .5*wi(jj,1)
        xi(jj,n+1) = 0.
        xr(jj,n+1) = 0.
 204  continue
      do 105 i=3,n,2
cdir$ ivdep
      do 105 j=1,m
         jj=1+(j-1)*jmp
         xi(jj,i-1) = wr(jj,i+1)
         xi(jj,i) = xi(jj,i-2)-wr(jj,i)
         xr(jj,i-1) = -wi(jj,i+1)
         xr(jj,i) = xr(jj,i-2)+wi(jj,i)
  105 continue
      return
      end if
      end

cvd$g nolstval
      subroutine vchbb(x,w,n,m,inc,jmp,wsave)
      implicit none
      real wsave(2*n+15),x(inc,1),w(inc,1)
      integer n,m,inc,jmp,nm1,np1,ns2,j,jj,k,kc,i
      real x1p3,tx2,t1,t2
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n.gt.1) then
         if (n.eq.3) then
            do j=1,m
               jj=1+(j-1)*jmp
               x1p3 = x(jj,1)+x(jj,3)
               tx2 = x(jj,2)
               x(jj,2) = (x(jj,1)-x(jj,3))
               x(jj,1) = (x1p3+tx2)
               x(jj,3) = (x1p3-tx2)
            end do
         else
cdir$ ivdep
cvd$ nodepchk
            do j=1,m
               jj=1+(j-1)*jmp
               w(jj,1) = x(jj,1)-x(jj,n)
               x(jj,1) = x(jj,1)+x(jj,n)
            end do
            do k=2,ns2
               kc = np1-k
cdir$ ivdep
cvd$ nodepchk
               do j=1,m
                  jj=1+(j-1)*jmp
                  t1 = .5*(x(jj,k)+x(jj,kc))
                  t2 = x(jj,k)-x(jj,kc)
                  w(jj,1) = w(jj,1)+.5*wsave(kc)*t2
                  t2 = .5*wsave(k)*t2
                  x(jj,k) = t1-t2
                  x(jj,kc) = t1+t2
               end do
            end do
            do j=1,m
               jj=1+(j-1)*jmp
C               x(jj,ns2+1) = 2.*x(jj,ns2+1)
            end do
            call rfftf1 (x,x(1,2),w,w(1,2),nm1,m,
     &           2*inc,jmp,wsave(n+1),wsave(n+n))
            do j=1,m
               jj=1+(j-1)*jmp
               x(jj,2) = w(jj,1)
            end do
            do i=4,n,2
cdir$ ivdep
               do j=1,m
                  jj=1+(j-1)*jmp
                  x(jj,i) = x(jj,i-2)-x(jj,i)
               end do
            end do
         end if
      end if

      end subroutine vchbb

cvd$g nolstval
      subroutine vchbf(x,w,n,m,inc,jmp,wsave)
      implicit none
      real wsave(2*n+15),x(inc,1),w(inc,1)
      integer nm1,np1,ns2,j,jj,n,m,inc,jmp,i,kc,k
      real x1p3,tx2,t1,t2
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n.gt.1) then
         if (n.eq.3) then
            do j=1,m
               jj=1+(j-1)*jmp
               x1p3 = x(jj,1)+x(jj,3)
               tx2 = x(jj,2)+x(jj,2)
               x(jj,2) = (x(jj,1)-x(jj,3))*.5
               x(jj,1) = (x1p3+tx2)*.25
               x(jj,3) = (x1p3-tx2)*.25
            end do
         else
cdir$ ivdep
            do j=1,m
               jj=1+(j-1)*jmp
               w(jj,1) = (x(jj,1)-x(jj,n))*.5
               x(jj,1) = (x(jj,1)+x(jj,n))*.5
            end do
            do k=2,ns2
               kc = np1-k
cdir$ ivdep
               do j=1,m
                  jj=1+(j-1)*jmp
                  t1 = .5*(x(jj,k)+x(jj,kc))
                  t2 = x(jj,k)-x(jj,kc)
                  w(jj,1) = w(jj,1)+.5*wsave(kc)*t2
                  t2 = .5*wsave(k)*t2
                  x(jj,k) = t1-t2
                  x(jj,kc) = t1+t2
               end do
            end do
            call rfftf1 (x,x(1,2),w,w(1,2),nm1,m,
     &           2*inc,jmp,wsave(n+1),wsave(n+n))
            do j=1,m
               jj=1+(j-1)*jmp
               x(jj,2) = w(jj,1)
            end do
            do i=4,n,2
               do j=1,m
                  jj=1+(j-1)*jmp
                  x(jj,i) = x(jj,i-2)-x(jj,i)
               end do
            end do
cdir$ ivdep
            do j=1,m
               jj=1+(j-1)*jmp
               x(jj,1) = x(jj,1)*.5
               x(jj,n) = x(jj,n)*.5
            end do
         end if
      end if

      end subroutine vchbf

cvd$g nolstval
      subroutine vcost (x,w,n,m,inc,jmp,wsave)
      implicit none
      real wsave(2*n+15),x(inc,1),w(inc,1)
      integer n,m,inc,jmp
      integer nm1,np1,ns2,j,jj,k,kc,i
      real x1p3,tx2,t1,t2
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n.gt.1) then
        if (n.eq.3) then
          do 101 j=1,m
            jj=1+(j-1)*jmp
            x1p3 = x(jj,1)+x(jj,3)
            tx2 = x(jj,2)+x(jj,2)
            x(jj,2) = x(jj,1)-x(jj,3)
            x(jj,1) = x1p3+tx2
            x(jj,3) = x1p3-tx2
 101        continue
        else
cdir$ ivdep
          do 102 j=1,m
            jj=1+(j-1)*jmp
            w(jj,1) = x(jj,1)-x(jj,n)
            x(jj,1) = x(jj,1)+x(jj,n)
 102      continue
          do 103 k=2,ns2
            kc = np1-k
cdir$ ivdep
            do 103 j=1,m
              jj=1+(j-1)*jmp
              t1 = x(jj,k)+x(jj,kc)
              t2 = x(jj,k)-x(jj,kc)
              w(jj,1) = w(jj,1)+wsave(kc)*t2
              x(jj,k) = t1-wsave(k)*t2
              x(jj,kc) = t1+wsave(k)*t2
 103      continue
          do 104 j=1,m
            jj=1+(j-1)*jmp
            x(jj,ns2+1) = 2.*x(jj,ns2+1)
 104      continue
          call rfftf1 (x,x(1,2),w,w(1,2),nm1,m,
     &      2*inc,jmp,wsave(n+1),wsave(n+n))
          do 105 j=1,m
            jj=1+(j-1)*jmp
c            w(jj,2) = x(jj,3)
            x(jj,2) = w(jj,1)
 105      continue
          do 106 i=4,n,2
cdir$ ivdep
          do 106 j=1,m
            jj=1+(j-1)*jmp
            x(jj,i) = x(jj,i-2)-x(jj,i)
 106      continue
          do 107 j=1,m
            jj=1+(j-1)*jmp
 107     continue
        end if
      end if
      return
      end

      subroutine vcosti (n,wsave,ifail)
      implicit none
      real wsave(2*n+15)
      real pi,dt,fk
      integer n,ifail,nm1,np1,ns2,kc,k,ifailn
      parameter (pi=3.14159265358979)

      ifailn=0
      if (n .eq. 2) ifailn=1
      if (n .gt. 3) then
        nm1 = n-1
        np1 = n+1
        ns2 = n/2
        dt = pi/real(nm1)
        fk = 0.
        do 101 k=2,ns2
          kc = np1-k
          fk = fk+1.
          wsave(k) = 2.*sin(fk*dt)
          wsave(kc) = 2.*cos(fk*dt)
  101   continue
        call rffti1 (nm1,wsave(n+1),wsave(n+n),ifailn)
      end if
      if (ifail.eq.0.and.ifailn.ne.0) then
          if (ifailn.eq.1) then
            write(*,*) 'VCOSTI-ERROR-N MUST BE ODD'
          else
            write(*,*) 'VCOSTI-ERROR-CANNOT FACTOR',ifailn
          end if
        stop
      end if
      if (ifail.ne.0) ifail=ifailn
      return
      end

cvd$g nolstval
      subroutine vrfftb (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n
      real cr(1+inc*(n/2)+jmp*(m-1)),ci(1+inc*(n/2)+jmp*(m-1))
      real wr(1+inc*(n/2)+jmp*(m-1)),wi(1+inc*(n/2)+jmp*(m-1))
      real wsave(n+15)
      if (n .eq. 1) return
      call rfftb1 (cr,ci,wr,wi,n,m,inc,jmp,wsave,wsave(n+1))
      return
      end



      subroutine rfftb1 (cr,ci,wr,wi,n,m,inc,jmp,wa,ifac)
      implicit none
      integer inc,jmp,n,m
      integer nf,j,jj,i,na,l1,iw,k1,ip,l2,ido,idl1,ix2,ix3,ix4
      real cr(1+inc*(n/2)+jmp*(m-1)),ci(1+inc*(n/2)+jmp*(m-1))
      real wr(1+inc*(n/2)+jmp*(m-1)),wi(1+inc*(n/2)+jmp*(m-1))
      real cl(8196)
      integer ifac(15)
      real wa(n)
      nf = ifac(2)
      do 1111 j=1,m
        jj=1+jmp*(j-1)
      do 1112 i=1,n/2+1
      cl(i)=cr(jj+(i-1)*inc)
      cl(n/2+1+i)=ci(jj+(i-1)*inc)
1112  continue
      na = 0
      l1 = 1
      iw = 1
        cl(n/2+2)=cl(1)
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (ido.eq.1) then
           if (na.eq.0) then
             if (ip.eq.4)call radb4l (l1,inc,
     &cl,cl(n/2+2),cr(jj),ci(jj))
             if (ip.eq.2)call radb2l (l1,inc,
     &cl,cl(n/2+2),cr(jj),ci(jj))
           else
             if (ip.eq.4)call radb4l (l1,inc,
     &cl(n+3),cl(n+n/2+4),cr(jj),ci(jj))
             if (ip.eq.2)call radb2l (l1,inc,
     &cl(n+3),cl(n+n/2+4),cr(jj),ci(jj))
           end if
         else
           if (na .eq. 0) then
             if (ip.eq.4)call radb4 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw),wa(ix2),wa(ix3))
             if (ip.eq.2)call radb2 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw))
             if (ip.eq.3)call radb3 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw),wa(ix2))
             if (ip.eq.5)call radb5 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                wa(iw),wa(ix2),wa(ix3),wa(ix4))
           else
             if (ip.eq.4)call radb4 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &                   wa(iw),wa(ix2),wa(ix3))
             if (ip.eq.2)call radb2 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &                   wa(iw))
             if (ip.eq.3)call radb3 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &          wa(iw),wa(ix2))
             if (ip.eq.5)call radb5 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &          wa(iw),wa(ix2),wa(ix3),wa(ix4))
           end if
         end if
         if (ip.gt.5) then
           write(*,*) 'factor ',ip,' not implemented'
           stop
         end if
         na=1-na
         l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
cc      if (na .ne. 0) then
cc        do 117 i=1,n/2
cc        do 117 j=1,m
cc          jj=1+jmp*(j-1)
cc          cr(jj+inc*(i-1)) = cl(n+2+i)
cc          ci(jj+inc*(i-1)) = cl(n+n/2+3+i)
cc  117   continue
cc      else
cc        do 119 i=1,n/2
cc          cr(jj+inc*(i-1)) = cl(i)
cc          ci(jj+inc*(i-1)) = cl(n/2+1+i)
cc  119   continue
cc      end if
cc      do 118 j=1,m
cc        jj=1+jmp*(j-1)
        cr(jj+inc*(n/2)) = 0.0
        ci(jj+inc*(n/2)) = 0.0
cc 118  continue
1111  continue
      return
      end

      subroutine radb2l (l1,inc,ccr,cci,chr,chi)
      implicit none
      integer k,l1,inc
      real ccr(l1+1),cci(l1),chr(inc,l1+1),chi(inc,l1)
      if (mod(l1,2).eq.0) then
        do 101 k=1,l1/2
cc        do 101 j=1,m
cc          jj=1+jmp*(j-1)
          chr(1,k) = cci(2*k-1)+ccr(2*k)
          chr(1,k+l1/2) = cci(2*k-1)-ccr(2*k)
          chi(1,k) = cci(2*k)+ccr(2*k+1)
          chi(1,k+l1/2) = cci(2*k)-ccr(2*k+1)
 101    continue
      else
        do 102 k=1,l1/2
cc        do 102 j=1,m
cc          jj=1+jmp*(j-1)
          chr(1,k) = cci(2*k-1)+ccr(2*k)
          chi(1,k+l1/2) = cci(2*k-1)-ccr(2*k)
          chi(1,k) = cci(2*k)+ccr(2*k+1)
          chr(1,k+l1/2+1) = cci(2*k)-ccr(2*k+1)
 102      continue
cc        do 103 j=1,m
cc          jj=1+jmp*(j-1)
          chr(1,l1/2+1) = cci(l1)+ccr(l1+1)
          chi(1,l1) = cci(l1)-ccr(l1+1)
cc 103    continue
      end if
      return
      end

      subroutine radb2 (ido,l1,ccr,cci,chr,chi,wa1)
      implicit none
      real ccr(ido/2,2,l1),cci(ido/2,2,l1)
      real chr(ido/2,l1,2),chi(ido/2,l1,2)
      real wa1(1)
      integer ic,k,ido,l1,idp2,i
      real tr2,ti2
      do 101 k=1,l1
cc      do 101 j=1,m
cc        jj=1+jmp*(j-1)
        chi(1,k,1) = cci(1,1,k)+ccr(ido/2+1,2,k)
        chi(1,k,2) = cci(1,1,k)-ccr(ido/2+1,2,k)
  101 continue
      if (ido.gt.2) then
        idp2 = ido/2+2
        do 103 k=1,l1
        do 103 i=2,ido/2
          ic = idp2-i
cc          do 104 j=1,m
cc            jj=1+jmp*(j-1)
            chr(i,k,1) = ccr(i,1,k)+ccr(ic,2,k)
            tr2 = ccr(i,1,k)-ccr(ic,2,k)
            chi(i,k,1) = cci(i,1,k)-cci(ic,2,k)
            ti2 = cci(i,1,k)+cci(ic,2,k)
            chr(i,k,2) = wa1(2*i-3)*tr2-wa1(2*i-2)*ti2
            chi(i,k,2) = wa1(2*i-3)*ti2+wa1(2*i-2)*tr2
cc  104     continue
  103   continue
      end if
      do 106 k=1,l1
cc      do 106 j=1,m
cc        jj=1+jmp*(j-1)
        chr(ido/2+1,k,1) = ccr(ido/2+1,1,k)+ccr(ido/2+1,1,k)
        chr(ido/2+1,k,2) = -(cci(1,2,k)+cci(1,2,k))
  106 continue
      return
      end

      subroutine radb3 (ido,l1,ccr,cci,chr,chi,wa1,wa2)
      implicit none
      real ccr(ido/2,3,l1),cci(ido/2,3,l1)
      real chr(ido/2,l1,3),chi(ido/2,l1,3)
      real  wa1(1),wa2(1)
      real taur,taui
      parameter(taur=-.5,taui=.866025403784439)
      integer k,l1,ido,ic,i,idp2
      real tr2,cr2,ci3,dr2,di2,ci2
      real tr1,ti2,cr3
      do 101 k=1,l1
cc      do 101 j=1,m
cc         jj=1+jmp*(j-1)
         tr2 = ccr(ido/2+1,2,k)+ccr(ido/2+1,2,k)
         cr2 = cci(1,1,k)+taur*tr2
         chi(1,k,1) = cci(1,1,k)+tr2
         ci3 = taui*(cci(1,3,k)+cci(1,3,k))
         chi(1,k,2) = cr2-ci3
         chi(1,k,3) = cr2+ci3
  101 continue
      idp2 = ido/2+2
      do 103 k=1,l1
      do 103 i=2,ido/2
        ic = idp2-i
cc        do 104 j=1,m
cc          jj=1+jmp*(j-1)
          tr2 = ccr(i,3,k)+ccr(ic,2,k)
          cr3 = taui*(ccr(i,3,k)-ccr(ic,2,k))
          cr2 = ccr(i,1,k)+taur*tr2
          chr(i,k,1) = ccr(i,1,k)+tr2
          ti2 = cci(i,3,k)-cci(ic,2,k)
          ci3 = taui*(cci(i,3,k)+cci(ic,2,k))
          ci2 = cci(i,1,k)+taur*ti2
          chi(i,k,1) = cci(i,1,k)+ti2
          dr2 = cr2-ci3
          di2 = ci2+cr3
          chr(i,k,2) = wa1(2*i-3)*dr2-wa1(2*i-2)*di2
          chi(i,k,2) = wa1(2*i-3)*di2+wa1(2*i-2)*dr2
C          di3 = ci2-cr3
C          dr3 = cr2+ci3
          chr(i,k,3) = wa2(2*i-3)*(cr2+ci3)-wa2(2*i-2)*(ci2-cr3)
          chi(i,k,3) = wa2(2*i-3)*(ci2-cr3)+wa2(2*i-2)*(cr2+ci3)
cc  104   continue
  103 continue
      do 106 k=1,l1
cc      do 106 j=1,m
cc        jj=1+jmp*(j-1)
        tr1 = ccr(ido/2+1,1,k)-ccr(ido/2+1,3,k)
        chr(ido/2+1,k,1) = 2.*ccr(ido/2+1,1,k)+ccr(ido/2+1,3,k)
        chr(ido/2+1,k,2) = tr1-(2.*taui)*cci(1,2,k)
        chr(ido/2+1,k,3) = -tr1-(2.*taui)*cci(1,2,k)
  106 continue
      return
      end

      subroutine radb4l (l1,inc,ccr,cci,chr,chi)
      implicit none
      real ccr(2,l1+1),cci(2,l1)
      real chr(inc,(l1+1)/2,4),chi(inc,(l1+1)/2,4)
      integer l1,k,inc
      real tr1,tr2,tr3,tr4
      if (mod(l1,2).eq.0) then
        do 101 k=1,l1/2
cc        do 101 j=1,m
cc           jj=1+jmp*(j-1)
C          tr1 = cci(1,2*k-1)-ccr(1,2*k)
C          tr2 = cci(1,2*k-1)+ccr(1,2*k)
C          tr3 = 2.*ccr(2,2*k-1)
C          tr4 = 2.*cci(2,2*k-1)
          chr(1,k,4) = cci(1,2*k-1)-ccr(1,2*k)+2.*cci(2,2*k-1)
          chr(1,k,2) = cci(1,2*k-1)-ccr(1,2*k)-2.*cci(2,2*k-1)
          chr(1,k,1) = cci(1,2*k-1)+ccr(1,2*k)+2.*ccr(2,2*k-1)
          chr(1,k,3) = cci(1,2*k-1)+ccr(1,2*k)-2.*ccr(2,2*k-1)
C          ti1 = cci(1,2*k)-ccr(1,2*k+1)
C          ti2 = cci(1,2*k)+ccr(1,2*k+1)
C          ti3 = 2.*ccr(2,2*k)
C          ti4 = 2.*cci(2,2*k)
          chi(1,k,4) = cci(1,2*k)-ccr(1,2*k+1)+2.*cci(2,2*k)
          chi(1,k,2) = cci(1,2*k)-ccr(1,2*k+1)-2.*cci(2,2*k)
          chi(1,k,1) = cci(1,2*k)+ccr(1,2*k+1)+2.*ccr(2,2*k)
          chi(1,k,3) = cci(1,2*k)+ccr(1,2*k+1)-2.*ccr(2,2*k)
 101    continue
      else
        do 102 k=1,l1/2
cc        do 102 j=1,m
cc          jj=1+jmp*(j-1)
          tr1 = cci(1,2*k-1)-ccr(1,2*k)
          tr2 = cci(1,2*k-1)+ccr(1,2*k)
          tr3 = ccr(2,2*k-1)+ccr(2,2*k-1)
          tr4 = cci(2,2*k-1)+cci(2,2*k-1)
          chr(1,k,1) = tr2+tr3
          chr(1,k-1,3) = tr2-tr3
          chi(1,k-2,4) = tr1+tr4
          chi(1,k-1,2) = tr1-tr4
          tr1 = cci(1,2*k)-ccr(1,2*k+1)
          tr2 = cci(1,2*k)+ccr(1,2*k+1)
          tr3 = ccr(2,2*k)+ccr(2,2*k)
          tr4 = cci(2,2*k)+cci(2,2*k)
          chi(1,k,1) = tr2+tr3
          chi(1,k-1,3) = tr2-tr3
          chr(1,k-1,4) = tr1+tr4
          chr(1,k,2) = tr1-tr4
 102    continue
cc        do 103 j=1,m
cc          jj=1+jmp*(j-1)
          tr1 = cci(1,l1)-ccr(1,l1+1)
          tr2 = cci(1,l1)+ccr(1,l1+1)
          tr3 = 2.*ccr(2,l1)
          tr4 = 2.*cci(2,l1)
          chr(1,l1/2+1,1) = tr2+tr3
          chr(1,l1/2,3) = tr2-tr3
          chi(1,l1/2-1,4) = tr1+tr4
          chi(1,l1/2,2) = tr1-tr4
cc 103    continue
      end if
      return
      end

      subroutine radb4 (ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3)
      implicit none
      real ccr(ido/2,4,l1),cci(ido/2,4,l1)
      real chr(ido/2,l1,4),chi(ido/2,l1,4)
      real wa1(1),wa2(1),wa3(1)
      real sqrt2
      parameter (sqrt2=1.414213562373095)
      integer k,l1,idp2,i,ic,ido
      real tr1,tr2,tr3,ti4,cr3,ti1,ti2,tr4,ti3
      do 101 k=1,l1
cc      do 101 j=1,m
cc        jj=1+jmp*(j-1)
        tr1 = cci(1,1,k)-ccr(ido/2+1,4,k)
        tr2 = cci(1,1,k)+ccr(ido/2+1,4,k)
C        tr3 = 2.*ccr(ido/2+1,2,k)
C        tr4 = 2.*cci(1,3,k)
        chi(1,k,4) = tr1+2.*cci(1,3,k)
        chi(1,k,2) = tr1-2.*cci(1,3,k)
        chi(1,k,1) = tr2+2.*ccr(ido/2+1,2,k)
        chi(1,k,3) = tr2-2.*ccr(ido/2+1,2,k)
  101 continue
      if (ido.gt.2) then
        idp2 = ido/2+2
        do 104 k=1,l1
        do 104 i=2,ido/2
          ic = idp2-i
cc          do 103 j=1,m
cc            jj=1+jmp*(j-1)
            tr3 = ccr(i,3,k)+ccr(ic,2,k)
            ti4 = ccr(i,3,k)-ccr(ic,2,k)
            tr2 = ccr(i,1,k)+ccr(ic,4,k)
            cr3 = tr2-tr3
            tr1 = ccr(i,1,k)-ccr(ic,4,k)
            chr(i,k,1) = tr2+tr3
            ti1 = cci(i,1,k)+cci(ic,4,k)
            ti2 = cci(i,1,k)-cci(ic,4,k)
            ti3 = cci(i,3,k)-cci(ic,2,k)
            chi(i,k,1) = ti2+ti3
            tr4 = cci(i,3,k)+cci(ic,2,k)
C            ci3 = ti2-ti3
            chr(i,k,3) = wa2(2*i-3)*cr3-wa2(2*i-2)*(ti2-ti3)
            chi(i,k,3) = wa2(2*i-3)*(ti2-ti3)+wa2(2*i-2)*cr3
C            cr2 = tr1-tr4
C            cr4 = tr1+tr4
C            ci2 = ti1+ti4
C            ci4 = ti1-ti4
            chr(i,k,2) = wa1(2*i-3)*(tr1-tr4)-wa1(2*i-2)*(ti1+ti4)
            chi(i,k,2) = wa1(2*i-3)*(ti1+ti4)+wa1(2*i-2)*(tr1-tr4)
            chr(i,k,4) = wa3(2*i-3)*(tr1+tr4)-wa3(2*i-2)*(ti1-ti4)
            chi(i,k,4) = wa3(2*i-3)*(ti1-ti4)+wa3(2*i-2)*(tr1+tr4)
cc  103     continue
  104   continue
      end if
      do 106 k=1,l1
cc      do 106 j=1,m
cc         jj=1+jmp*(j-1)
         ti1 = cci(1,2,k)+cci(1,4,k)
C         ti2 = cci(1,4,k)-cci(1,2,k)
         chr(ido/2+1,k,3) = 2.*(cci(1,4,k)-cci(1,2,k))
         tr1 = ccr(ido/2+1,1,k)-ccr(ido/2+1,3,k)
C         tr2 = ccr(ido/2+1,1,k)+ccr(ido/2+1,3,k)
         chr(ido/2+1,k,1) =
     &        2.*(ccr(ido/2+1,1,k)+ccr(ido/2+1,3,k))
         chr(ido/2+1,k,2) = sqrt2*(tr1-ti1)
         chr(ido/2+1,k,4) = -sqrt2*(tr1+ti1)
  106 continue
      return
      end


      subroutine radb5 (ido,l1,ccr,cci,chr,chi,
     &     wa1,wa2,wa3,wa4)
      implicit none
      real ccr(ido/2,5,l1),cci(ido/2,5,l1)
      real chr(ido/2,l1,5),chi(ido/2,l1,5)
      real wa1(1)     ,wa2(1)     ,wa3(1)     ,wa4(1)
      real tr11,ti11,tr12,ti12,tr13
      integer k,l1,ido,i,idp2,ic
      real ti5,ti4,tr2,tr3,cr2,cr3,ci5
      real t1,t2,t3,t4,t5,t6,t7
      real ti2,ti3,ci2,ci3,tr5,tr4
      real ci4,cr5,cr4,dr3,dr4,di3,di4,dr5,dr2,di5,di2
      parameter(tr11= .309016994374947)
      parameter(ti11= .951056516295155)
      parameter(tr12=-.809016994374947)
      parameter(ti12= .587785252292473)
      parameter(tr13= .55901699437495)

      do 101 k=1,l1
cc      do 101 j=1,m
cc         jj=1+jmp*(j-1)
         ti5 = cci(1,3,k)+cci(1,3,k)
         ti4 = cci(1,5,k)+cci(1,5,k)
         tr2 = ccr(ido/2+1,2,k)+ccr(ido/2+1,2,k)
         tr3 = ccr(ido/2+1,4,k)+ccr(ido/2+1,4,k)
         chi(1,k,1) = cci(1,1,k)+tr2+tr3
         cr2 = cci(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cci(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         chi(1,k,2) = cr2-ci5
         chi(1,k,3) = cr3-ci4
         chi(1,k,4) = cr3+ci4
         chi(1,k,5) = cr2+ci5
  101 continue
      idp2 = ido/2+2
      do 103 k=1,l1
      do 103 i=2,ido/2
        ic = idp2-i
cc        do 102 j=1,m
cc            jj=1+jmp*(j-1)
            ti5 = cci(i,3,k)+cci(ic,2,k)
            ti2 = cci(i,3,k)-cci(ic,2,k)
            ti4 = cci(i,5,k)+cci(ic,4,k)
            ci5 = ti11*ti5+ti12*ti4
            ci4 = ti12*ti5-ti11*ti4
            ti3 = cci(i,5,k)-cci(ic,4,k)
            ci2 = cci(i,1,k)+tr11*ti2+tr12*ti3
            ci3 = cci(i,1,k)+tr12*ti2+tr11*ti3
            chi(i,k,1) = cci(i,1,k)+ti2+ti3
            tr5 = ccr(i,3,k)-ccr(ic,2,k)
            tr2 = ccr(i,3,k)+ccr(ic,2,k)
            tr4 = ccr(i,5,k)-ccr(ic,4,k)
            cr5 = ti11*tr5+ti12*tr4
            cr4 = ti12*tr5-ti11*tr4
            tr3 = ccr(i,5,k)+ccr(ic,4,k)
            cr2 = ccr(i,1,k)+tr11*tr2+tr12*tr3
            cr3 = ccr(i,1,k)+tr12*tr2+tr11*tr3
            chr(i,k,1) = ccr(i,1,k)+tr2+tr3
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            chr(i,k,3) = wa2(2*i-3)*dr3-wa2(2*i-2)*di3
            chi(i,k,3) = wa2(2*i-3)*di3+wa2(2*i-2)*dr3
            chr(i,k,4) = wa3(2*i-3)*dr4-wa3(2*i-2)*di4
            chi(i,k,4) = wa3(2*i-3)*di4+wa3(2*i-2)*dr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            chr(i,k,2) = wa1(2*i-3)*dr2-wa1(2*i-2)*di2
            chi(i,k,2) = wa1(2*i-3)*di2+wa1(2*i-2)*dr2
            chr(i,k,5) = wa4(2*i-3)*dr5-wa4(2*i-2)*di5
            chi(i,k,5) = wa4(2*i-3)*di5+wa4(2*i-2)*dr5
cc  102     continue
  103 continue
      do 106 k=1,l1
cc      do 106 j=1,m
cc         jj=1+jmp*(j-1)
         t1 = ccr(ido/2+1,1,k)+ccr(ido/2+1,3,k)
         t2 = .5*t1-ccr(ido/2+1,5,k)
         t3 = (2.*tr13)*(ccr(ido/2+1,1,k)-ccr(ido/2+1,3,k))
         t4 = (2.*ti12)*cci(1,2,k)+(2.*ti11)*cci(1,4,k)
         t5 = (2.*ti11)*cci(1,2,k)-(2.*ti12)*cci(1,4,k)
         t6 = t3+t2
         t7 = t3-t2
         chr(ido/2+1,k,1) = ccr(ido/2+1,5,k)+2.*t1
         chr(ido/2+1,k,2) = t6-t4
         chr(ido/2+1,k,3) = t7-t5
         chr(ido/2+1,k,4) = -t7-t5
         chr(ido/2+1,k,5) = -t6-t4
  106 continue
      return
      end

cvd$g nolstval
      subroutine vrfftf (cr,ci,wr,wi,n,m,inc,jmp,wsave)
      implicit none
      integer inc,jmp,m,n
      real cr(1+inc*(n/2)+jmp*(m-1)),ci(1+inc*(n/2)+jmp*(m-1))
      real wr(1+inc*(n/2)+jmp*(m-1)),wi(1+inc*(n/2)+jmp*(m-1))
      real wsave(n+15)
      integer j
      if (n .eq. 1) return
      call rfftf1 (cr,ci,wr,wi,n,m,inc,jmp,wsave,wsave(n+1))
      do 117 j=1,m
        ci(1+jmp*(j-1)+inc*n/2)=0.
 117  continue
      return
      end


      subroutine rfftf1 (cr,ci,wr,wi,n,m,inc,jmp,wa,ifac)
      implicit none
      integer n,m,inc,jmp,ifac(15)
      real wa(n)
      real cl(8196)
      real cr(1+inc*(n/2)+jmp*(m-1)),ci(1+inc*(n/2)+jmp*(m-1))
      real wr(1+inc*(n/2)+jmp*(m-1)),wi(1+inc*(n/2)+jmp*(m-1))
      integer nf,na,l2,iw,k1,kh,ido,idl1,ip,l1,ix2,ix3,ix4,i,j,ii
      integer jj
      nf = ifac(2)
      do 1111 j=1,m
        jj=1+jmp*(j-1)
c      do 1112 i=1,n/2
c      cl(i)=cr(jj+(i-1)*inc)
c      cl(i+n/2+1)=ci(jj+(i-1)*inc)
c1112  continue
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (ido.eq.1) then
           if (ip.eq.4)call radf4f (l1,inc,
     &cr(jj),ci(jj),cl(n+3),cl(n+n/2+4))
           if (ip.eq.2)call radf2f (l1,inc,
     &cr(jj),ci(jj),cl(n+3),cl(n+n/2+4))
         else
           if (na .eq. 0) then
             if (ip.eq.4)call radf4 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw),wa(ix2),wa(ix3))
             if (ip.eq.2)call radf2 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw))
             if (ip.eq.3)call radf3 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                   wa(iw),wa(ix2))
             if (ip.eq.5)call radf5 (ido,l1,
     &cl,cl(n/2+2),cl(n+3),cl(n+n/2+4),
     &                wa(iw),wa(ix2),wa(ix3),wa(ix4))
           else
             if (ip.eq.4)call radf4 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &                   wa(iw),wa(ix2),wa(ix3))
             if (ip.eq.2)call radf2 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &                   wa(iw))
             if (ip.eq.3)call radf3 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &          wa(iw),wa(ix2))
             if (ip.eq.5)call radf5 (ido,l1,
     &cl(n+3),cl(n+n/2+4),cl,cl(n/2+2),
     &          wa(iw),wa(ix2),wa(ix3),wa(ix4))
           end if
         end if
           if (ip.gt.5) then
             write(*,*) 'Factor ',ip,' not implemented'
             stop
           end if
        l2 = l1
  111 continue
      if (na .eq. 0) then
        do 112 i=2,n/2
cc        do 112 j=1,m
          ii=jj+inc*(i-1)
          cr(ii) = cl(n+2+i)
          ci(ii) = cl(n+n/2+3+i)
  112   continue
cc        do 116 j=1,m
          cr(jj)=cl(n+n/2+4)
cc 116   continue
cc        do 113 j=1,m
          cr(jj+inc*n/2)=cl(n+n/2+3)
cc 113    continue
      else
        do 212 i=2,n/2
cc        do 112 j=1,m
          ii=jj+inc*(i-1)
          cr(ii) = cl(i)
          ci(ii) = cl(n/2+1+i)
  212   continue
cc        do 116 j=1,m
          cr(jj)=cl(n/2+2)
cc 116   continue
cc        do 113 j=1,m
          cr(jj+inc*n/2)=cl(n/2+1)
cc 113    continue
cc        do 114 j=1,m
cc          cr(jj)=ci(jj)
cc 114      continue
      end if
cc        do 115 j=1,m
          ci(jj)=0.
cc 115    continue
1111  continue
      return
      end


      subroutine radf2f (l1,inc,ccr,cci,chr,chi)
      implicit none
      integer inc,l1,k
      real ccr(inc,1),cci(inc,1)
      real  chr(2,1),chi(2,1)
      if (mod(l1,2).eq.0) then
      do 101 k=1,l1/2
cc      do 101 j=1,m
cc        jj=1+jmp*(j-1)
        chi(1,k) = ccr(1,k)+ccr(1,k+l1/2)
        chr(2,k) = ccr(1,k)-ccr(1,k+l1/2)
        chi(2,k) = cci(1,k)+cci(1,k+l1/2)
        chr(1,k+1) = cci(1,k)-cci(1,k+l1/2)
 101  continue
      else
      do 102 k=1,l1/2
cc      do 102 j=1,m
cc        jj=1+jmp*(j-1)
        chi(1,k) = ccr(1,k)+cci(1,k+l1/2)
        chr(2,k) = ccr(1,k)-cci(1,k+l1/2)
        chi(2,k) = cci(1,k)+ccr(1,k+l1/2+1)
        chr(1,k+1) = cci(1,k)-ccr(1,k+l1/2+1)
 102    continue
cc      do 103 j=1,m
cc        jj=1+jmp*(j-1)
        chi(1,l1/2+1) = ccr(1,l1/2+1)+cci(1,l1)
        chr(2,l1/2+1) = ccr(1,l1/2+1)-cci(1,l1)
cc 103  continue
      end if
      return
      end

      subroutine radf2 (ido,l1,ccr,cci,chr,chi,wa1)
      implicit none
      integer k,idp2,i,ic,ido,l1
      real tr2,ti2
      real chr(ido/2,2,l1),chi(ido/2,2,l1)
      real wa1(ido-1),ccr(ido/2,l1,2),cci(ido/2,l1,2)
      do 101 k=1,l1
cc      do 101 j=1,m
cc         jj=1+jmp*(j-1)
         chi(1,1,k) = cci(1,k,1)+cci(1,k,2)
         chr(ido/2+1,2,k) = cci(1,k,1)-cci(1,k,2)
  101 continue
      if (ido.gt.2) then
        idp2 = ido/2+2
        do 104 k=1,l1
        do 104 i=2,ido/2
          ic = idp2-i
cc          do 105 j=1,m
cc            jj=1+jmp*(j-1)
            tr2 = wa1(2*i-3)*ccr(i,k,2)+wa1(2*i-2)*cci(i,k,2)
            ti2 = wa1(2*i-3)*cci(i,k,2)-wa1(2*i-2)*ccr(i,k,2)
            chi(i,1,k) = cci(i,k,1)+ti2
            chi(ic,2,k) = ti2-cci(i,k,1)
            chr(i,1,k) = ccr(i,k,1)+tr2
            chr(ic,2,k) = ccr(i,k,1)-tr2
cc 105     continue
 104  continue
      end if
      do 106 k=1,l1
cc      do 106 j=1,m
cc        jj=1+jmp*(j-1)
        chi(1,2,k) = -ccr(ido/2+1,k,2)
        chr(ido/2+1,1,k) = ccr(ido/2+1,k,1)
106   continue
      return
      end

      subroutine radf3 (ido,l1,ccr,cci,chr,chi,wa1,wa2)
      implicit none
      integer ido,l1,k,i,ic,idp2
      real t1,ti3,tr3,tr2,ti2,di3,dr3,di2,dr2,cr2
      real  chr(ido/2,3,l1),chi(ido/2,3,l1)
      real  ccr(ido/2,l1,3),cci(ido/2,l1,3), wa1(1),wa2(1)
      real taur,taui

      parameter(taur=-.5,taui=.866025403784439)
      do 101 k=1,l1
cc      do 101 j=1,m
cc          jj=1+jmp*(j-1)
          cr2 = cci(1,k,2)+cci(1,k,3)
          chi(1,1,k) = cci(1,k,1)+cr2
          chi(1,3,k) = taui*(cci(1,k,3)-cci(1,k,2))
          chr(ido/2+1,2,k) = cci(1,k,1)+taur*cr2
  101 continue
      idp2 = ido/2+2
      do 102 k=1,l1
      do 102 i=2,ido/2
        ic = idp2-i
cc        do 103 j=1,m
cc          jj=1+jmp*(j-1)
          dr2 = wa1(2*i-3)*ccr(i,k,2)+wa1(2*i-2)*cci(i,k,2)
          di2 = wa1(2*i-3)*cci(i,k,2)-wa1(2*i-2)*ccr(i,k,2)
          dr3 = wa2(2*i-3)*ccr(i,k,3)+wa2(2*i-2)*cci(i,k,3)
          di3 = wa2(2*i-3)*cci(i,k,3)-wa2(2*i-2)*ccr(i,k,3)
C          cr2 = dr2+dr3
C          ci2 = di2+di3
          chr(i,1,k) = ccr(i,k,1)+(dr2+dr3)
          tr2 = ccr(i,k,1)+taur*(dr2+dr3)
          chi(i,1,k) = cci(i,k,1)+(di2+di3)
          ti2 = cci(i,k,1)+taur*(di2+di3)
          tr3 = taui*(di2-di3)
          ti3 = taui*(dr3-dr2)
          chr(i,3,k) = tr2+tr3
          chr(ic,2,k) = tr2-tr3
          chi(i,3,k) = ti2+ti3
          chi(ic,2,k) = ti3-ti2
cc 103    continue
 102  continue
      do 106 k=1,l1
cc      do 106 j=1,m
cc         jj=1+jmp*(j-1)
         t1 = ccr(ido/2+1,k,2)-ccr(ido/2+1,k,3)
         chr(ido/2+1,1,k) = .5*t1+ccr(ido/2+1,k,1)
         chr(ido/2+1,3,k) = ccr(ido/2+1,k,1)-t1
         chi(1,2,k) = -taui*(ccr(ido/2+1,k,2)+ccr(ido/2+1,k,3))
  106 continue
      return
      end



      subroutine radf4f (l1,inc,ccr,cci,chr,chi)
      implicit none
      integer l1,inc,k
      real tr1,tr2
      real ccr(inc,2*l1+1),cci(inc,2*l1),chr(4,l1+1),chi(4,l1)
      if (mod(l1,2).eq.0) then
        do 101 k=1,l1/2
cc        do 101 j=1,m
cc          jj=1+jmp*(j-1)
          tr1 = ccr(1,k+l1/2)+ccr(1,k+3*l1/2)
          chi(2,k) = ccr(1,k+3*l1/2)-ccr(1,k+l1/2)
          tr2 = ccr(1,k)+ccr(1,k+l1)
          chr(2,k) = ccr(1,k)-ccr(1,k+l1)
          chi(1,k) = tr1+tr2
          chr(3,k) = tr2-tr1
          tr1 = cci(1,k+l1/2)+cci(1,k+3*l1/2)
          chi(4,k) = cci(1,k+3*l1/2)-cci(1,k+l1/2)
          tr2 = cci(1,k)+cci(1,k+l1)
          chr(4,k) = cci(1,k)-cci(1,k+l1)
          chi(3,k) = tr1+tr2
          chr(1,k+1) = tr2-tr1
 101    continue
      else
        do 102 k=1,l1/2
cc        do 102 j=1,m
cc          jj=1+jmp*(j-1)
          tr1 = cci(1,k+l1/2)+cci(1,k+3*l1/2)
          chi(2,k) = cci(1,k+3*l1/2)-cci(1,k+l1/2)
          tr2 = ccr(1,k)+ccr(1,k+l1)
          chr(2,k) = ccr(1,k)-ccr(1,k+l1)
          chi(1,k) = tr1+tr2
          chr(3,k) = tr2-tr1
          tr1 = ccr(1,k+l1/2+1)+ccr(1,k+3*l1/2+1)
          chi(4,k) = ccr(1,k+3*l1/2+1)-ccr(1,k+l1/2+1)
          tr2 = cci(1,k)+cci(1,k+l1)
          chr(4,k) = cci(1,k)-cci(1,k+l1)
          chi(3,k) = tr1+tr2
          chr(1,k+1) = tr2-tr1
 102    continue
cc        do 103 j=1,m
cc          jj=1+jmp*(j-1)
          tr1 = cci(1,l1)+cci(1,2*l1)
          tr2 = ccr(1,l1/2+1)+ccr(1,3*l1/2+1)
          chi(1,l1/2+1) = tr1+tr2
          chr(3,l1/2+1) = tr2-tr1
          chr(2,l1/2+1) = ccr(1,l1/2+1)-ccr(1,3*l1/2+1)
          chi(2,l1/2+1) = cci(1,2*l1)-cci(1,l1)
cc 103    continue
      end if
      return
      end


      subroutine radf4 (ido,l1,ccr,cci,chr,chi,wa1,wa2,wa3)
      implicit none
      integer ido,l1,k,ic,i,idp2
      real tr1,tr2,ti1,ti3,ti2,cr3,ci3,ci4,cr4,ci2,cr2
      real ccr(ido/2,l1,4),cci(ido/2,l1,4)
      real chr(ido/2,4,l1),chi(ido/2,4,l1)
      real wa1(1)     ,wa2(1)     ,wa3(1)
      real hsqt2
      parameter (hsqt2=.7071067811865475)
      do 101 k=1,l1
cc      do 101 j=1,m
cc         jj=1+jmp*(j-1)
         tr1 = cci(1,k,2)+cci(1,k,4)
         tr2 = cci(1,k,1)+cci(1,k,3)
         chi(1,1,k) = tr1+tr2
         chr(ido/2+1,4,k) = tr2-tr1
         chr(ido/2+1,2,k) = cci(1,k,1)-cci(1,k,3)
         chi(1,3,k) = cci(1,k,4)-cci(1,k,2)
  101 continue
      if (ido.gt.2) then
  102   idp2 = ido/2+2
        do 104 k=1,l1
        do 104 i=2,ido/2
          ic = idp2-i
cc          do 103 j=1,m
cc            jj=1+jmp*(j-1)
            cr2 = wa1(2*i-3)*ccr(i,k,2)+wa1(2*i-2)*cci(i,k,2)
            ci2 = wa1(2*i-3)*cci(i,k,2)-wa1(2*i-2)*ccr(i,k,2)
            cr4 = wa3(2*i-3)*ccr(i,k,4)+wa3(2*i-2)*cci(i,k,4)
            ci4 = wa3(2*i-3)*cci(i,k,4)-wa3(2*i-2)*ccr(i,k,4)
C            tr1 = cr2+cr4
C            tr4 = cr4-cr2
C            ti1 = ci2+ci4
C            ti4 = ci2-ci4
            cr3 = wa2(2*i-3)*ccr(i,k,3)+wa2(2*i-2)*cci(i,k,3)
            ci3 = wa2(2*i-3)*cci(i,k,3)-wa2(2*i-2)*ccr(i,k,3)
            ti2 = cci(i,k,1)+ci3
            ti3 = cci(i,k,1)-ci3
            chi(i,1,k) = (ci2+ci4)+ti2
            chi(ic,4,k) = (ci2+ci4)-ti2
            chi(i,3,k) = (cr4-cr2)+ti3
            chi(ic,2,k) = (cr4-cr2)-ti3
C            tr2 = ccr(i,k,1)+cr3
C            tr3 = ccr(i,k,1)-cr3
            chr(i,1,k) = (cr2+cr4)+(ccr(i,k,1)+cr3)
            chr(ic,4,k) = (ccr(i,k,1)+cr3)-(cr2+cr4)
            chr(i,3,k) = (ci2-ci4)+(ccr(i,k,1)-cr3)
            chr(ic,2,k) = (ccr(i,k,1)-cr3)-(ci2-ci4)
cc  103     continue
  104   continue
      end if
      do 106 k=1,l1
cc      do 106 j=1,m
cc         jj=1+jmp*(j-1)
         ti1 = -hsqt2*(ccr(ido/2+1,k,2)+ccr(ido/2+1,k,4))
         tr1 = hsqt2*(ccr(ido/2+1,k,2)-ccr(ido/2+1,k,4))
         chr(ido/2+1,1,k) = tr1+ccr(ido/2+1,k,1)
         chr(ido/2+1,3,k) = ccr(ido/2+1,k,1)-tr1
         chi(1,2,k) = ti1-ccr(ido/2+1,k,3)
         chi(1,4,k) = ti1+ccr(ido/2+1,k,3)
  106 continue
  107 return
      end

      subroutine radf5 (ido,l1,ccr,cci,chr,chi,
     &     wa1,wa2,wa3,wa4)
      implicit none
      integer l1,ido,k,ic,i,idp2
      real dr2,di2,dr3,di3
      real cr3,ci4,ci5,cr2
      real ti5,tr4,tr5,ti3,ti2,ci3,t3,t4,t1,t2,ti4
      real t7,t6,t5,cr4,ci2,cr5,tr3,tr2,di5,dr5,di4,dr4
      real ccr(ido/2,l1,5),cci(ido/2,l1,5)
      real chr(ido/2,5,l1),chi(ido/2,5,l1)
      real wa1(1)     ,wa2(1)     ,wa3(1)     ,wa4(1)
      real tr11,ti11,tr12,ti12,tr13
      parameter(tr11= .309016994374947)
      parameter(ti11= .951056516295154)
      parameter(tr12=-.809016994374947)
      parameter(ti12= .587785252292473)
      parameter(tr13= .55901699437495)
      do 101 k=1,l1
cc      do 101 j=1,m
cc        jj=1+jmp*(j-1)
        cr2 = cci(1,k,5)+cci(1,k,2)
        ci5 = cci(1,k,5)-cci(1,k,2)
        ci4 = cci(1,k,4)-cci(1,k,3)
        chi(1,3,k) = ti11*ci5+ti12*ci4
        chi(1,5,k) = ti12*ci5-ti11*ci4
        cr3 = cci(1,k,4)+cci(1,k,3)
        chi(1,1,k) = cci(1,k,1)+cr2+cr3
        chr(ido/2+1,2,k) = cci(1,k,1)+tr11*cr2+tr12*cr3
        chr(ido/2+1,4,k) = cci(1,k,1)+tr12*cr2+tr11*cr3
 101  continue
      if (ido.gt.2) then
        idp2 = ido/2+2
        do 103 k=1,l1
        do 103 i=2,ido/2
          ic = idp2-i
cc          do 102 j=1,m
cc            jj=1+jmp*(j-1)
            dr2 = wa1(2*i-3)*ccr(i,k,2)+wa1(2*i-2)*cci(i,k,2)
            di2 = wa1(2*i-3)*cci(i,k,2)-wa1(2*i-2)*ccr(i,k,2)
            dr3 = wa2(2*i-3)*ccr(i,k,3)+wa2(2*i-2)*cci(i,k,3)
            di3 = wa2(2*i-3)*cci(i,k,3)-wa2(2*i-2)*ccr(i,k,3)
            dr4 = wa3(2*i-3)*ccr(i,k,4)+wa3(2*i-2)*cci(i,k,4)
            di4 = wa3(2*i-3)*cci(i,k,4)-wa3(2*i-2)*ccr(i,k,4)
            dr5 = wa4(2*i-3)*ccr(i,k,5)+wa4(2*i-2)*cci(i,k,5)
            di5 = wa4(2*i-3)*cci(i,k,5)-wa4(2*i-2)*ccr(i,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            chr(i,1,k) = ccr(i,k,1)+cr2+cr3
            tr2 = ccr(i,k,1)+tr11*cr2+tr12*cr3
            tr3 = ccr(i,k,1)+tr12*cr2+tr11*cr3
            cr5 = di2-di5
            ci2 = di2+di5
            cr4 = di3-di4
            ci3 = di3+di4
            chi(i,1,k) = cci(i,k,1)+ci2+ci3
            ti2 = cci(i,k,1)+tr11*ci2+tr12*ci3
            ti3 = cci(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            tr4 = ti12*cr5-ti11*cr4
            chr(i,5,k) = tr3+tr4
            chr(ic,4,k) = tr3-tr4
            chr(i,3,k) = tr2+tr5
            chr(ic,2,k) = tr2-tr5
            ti5 = ti11*ci5+ti12*ci4
            ti4 = ti12*ci5-ti11*ci4
            chi(i,3,k) = ti2+ti5
            chi(ic,2,k) = ti5-ti2
            chi(i,5,k) = ti3+ti4
            chi(ic,4,k) = ti4-ti3
  102     continue
  103   continue
      end if
      if (mod(ido,2) .eq. 1) return
      do 106 k=1,l1
cc      do 106 j=1,m
cc         jj=1+jmp*(j-1)
         t2 = ccr(ido/2+1,k,2)+ccr(ido/2+1,k,5)
         t1 = ccr(ido/2+1,k,2)-ccr(ido/2+1,k,5)
         t4 = ccr(ido/2+1,k,3)+ccr(ido/2+1,k,4)
         t3 = ccr(ido/2+1,k,3)-ccr(ido/2+1,k,4)
         t5 = t1-t3
         t6 = ccr(ido/2+1,k,1)+.25*t5
         t7 = tr13*(t1+t3)
         chr(ido/2+1,1,k) = t6+t7
         chr(ido/2+1,3,k) = t6-t7
         chr(ido/2+1,5,k) = ccr(ido/2+1,k,1)-t5
         chi(1,2,k) = -ti12*t2-ti11*t4
         chi(1,4,k) = -ti11*t2+ti12*t4
  106 continue
      return
      end

      subroutine vrffti (n,wsave,ifail)
      implicit none
      integer ifail,n,ifailn
      real wsave(n+15)
      if (n .eq. 1) then
        ifailn=1
      else
        call rffti1 (n,wsave(1),wsave(n+1),ifailn)
      end if
      if (ifail.eq.0.and.ifailn.ne.0) then
        if (ifailn.eq.1) then
          write(*,*) 'vrffti-error-n must be even'
        else
          write(*,*) 'vrffti-error-cannot factor',ifailn
        end if
        stop
      end if
      if (ifail.ne.0) ifail=ifailn
      return
      end


      subroutine rffti1 (n,wa,ifac,ifail)
      implicit none
      real wa(n)
      integer ifac(15),ntryh(4),nfac(4)
      data ntryh/2,3,4,5/
      real arg,fi,argld,tpi,argh
      integer is,nfm1,l1
      integer ii,ld,n,i,j,ipm,ido,l2,ip,k1,ifail,nl,nf,ntry,np
      parameter(tpi = 6.28318530717959)
      ifail=0
      nl = n
      nf = 0
c first factor n, select large factors first
      do 101 j=4,1,-1
        nfac(j)=0
        ntry=ntryh(j)
 102    if (nl.eq.nl/ntry*ntry) then
          nl=nl/ntry
          nf=nf+1
          nfac(j)=nfac(j)+1
          goto 102
        end if
 101  continue
      if (nl.ne.1) then
        ifail=nl
        return
      end if
c now insert the factors into ifac with the largest even factor last
      if (nfac(3).gt.0) then
        ifac(nf+2)=4
        nfac(3)=nfac(3)-1
      else
        if (nfac(1).gt.0) then
          ifac(nf+2)=2
          nfac(1)=nfac(1)-1
        else
          ifail=1
          return
        end if
      end if
      np=3
      do 103 j=1,4
        do 104 np=np,np-1+nfac(j)
          ifac(np)=ntryh(j)
 104    continue
 103  continue
      ifac(1) = n
      ifac(2) = nf
      argh = tpi/real(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = real(ld)*argh
            fi = 0.
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue
      return
      end

cvd$g nolstval
      subroutine vsint (x,w,n,m,inc,jmp,wsave)
      implicit none
      real x(1),wsave(1),w(1)
      integer np1,iw1,iw2,n,m,inc,jmp
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      call sint1(x,w,n,m,inc,jmp,wsave,wsave(iw1),wsave(iw2))
      return
      end

      subroutine sint1(x,w,n,m,inc,jmp,was,war,ifac)
      implicit none
      real war(1),was(1),x(inc,1),w(inc,1)
      integer ifac(1)
      integer n,m,jmp,j,jj,i,inc,kc,k,ns2,np1
      real t1,t2
      real sqrt3
      parameter(sqrt3=1.73205080756888)
      if (n.eq.1) then
        do 200 j=1,m
          jj=1+(j-1)*jmp
          x(jj,1) = x(jj,1)+x(jj,1)
 200    continue
      else
      np1 = n+1
      ns2 = n/2
      do 205 j=1,m
        jj=1+(j-1)*jmp
        w(jj,1) = 0.
 205  continue
      do 104 k=1,ns2
        kc = np1-k
cdir$ ivdep
        do 202 j=1,m
          jj=1+(j-1)*jmp
          t1 = x(jj,k)-x(jj,kc)
          t2 = was(k)*(x(jj,k)+x(jj,kc))
          w(jj,k+1) = t1+t2
          w(jj,kc+1) = t2-t1
  202   continue
  104 continue
      do 203 j=1,m
        jj=1+(j-1)*jmp
        w(jj,ns2+2) = 4.*x(jj,ns2+1)
 203  continue
      call rfftf1 (w,w(1,2),x,x(1,2),np1,m,2*inc,jmp,war,ifac)
      do 204 j=1,m
        jj=1+(j-1)*jmp
        x(jj,1) = .5*w(jj,1)
 204  continue
      do 105 i=3,n,2
cdir$ ivdep
      do 105 j=1,m
         jj=1+(j-1)*jmp
         x(jj,i-1) = -w(jj,i+1)
         x(jj,i) = x(jj,i-2)+w(jj,i)
  105 continue
      return
      end if
      end

      subroutine vsinti (n,wsave,ifail)
      implicit none
      real wsave(3*n/2+15)
      integer n,ifail,ifailn,ns2,np1,k
      real dt
      real pi
      parameter (pi=3.14159265358979)
      ifailn=0
      if (n .ne. 1) then
        ns2 = n/2
        np1 = n+1
        dt = pi/real(np1)
        do 101 k=1,ns2
          wsave(k) = 2.*sin(k*dt)
  101   continue
        call rffti1 (np1,wsave(ns2+1),wsave(ns2+n+2),ifailn)
        if (ifail.eq.0.and.ifailn.ne.0) then
          if (ifailn.eq.1) then
            write(*,*) 'VSINTI-ERROR-N MUST BE ODD'
          else
            write(*,*) 'VSINTI-ERROR-CANNOT FACTOR',ifailn
          end if
          stop
        end if
      end if
      if (ifail.ne.0) ifail=ifailn
      return
      end
