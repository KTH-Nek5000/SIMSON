c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wdisca(ur,ui,re,xl,zl,t,xs,dstar,fltype,
     &   bstart,blength,rlam,spanv,namnut,m,boxr,boxi,urx)
c
c     Writes m variables form ur to file namnut
c
      implicit none

      include 'par.f'

      character*80 namnut
      integer m,fltype
      real ur(memnx,memny,memnz,7),ui(memnx,memny,memnz,7)
      real re,xl,zl,t,xs,dstar
      real bstart,blength,rlam,spanv
      real urx(nx)
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real prey(2*ny+15)

      integer x,y,z,i,zb
c
c     Preprocess y-transform
c
      call vcosti(ny,prey,0)
c
c     Write file header
c
      open(unit=11,file=namnut,status='new',form='unformatted')
      write(11) re,.false.,xl,zl,t,xs
      write(11) nx,nyp,nzc,nfzsym
      write(11) fltype,dstar
      if (fltype.eq.-1) write(11) rlam
      if (fltype.eq.-2) write(11) rlam,spanv
      if (fltype.ge.4) write(11) bstart,blength,rlam,spanv

      do i=1,m
         do zb=1,nzc,mbz
            call getxy(boxr,boxi,zb,i,ur,ui)
            do z=zb,zb+mbz-1
               if (z+(i-1)*nzc.eq.nzc*m/100+1) write(*,*) 'Wrote 1%'
               if (z+(i-1)*nzc.eq.nzc*m/10+1) write(*,*) 'Wrote 10%'
               if (z+(i-1)*nzc.eq.nzc*m/2+1) write(*,*) 'Wrote 50%'
               if (z+(i-1)*nzc.eq.nzc*m*9/10+1) write(*,*) 'Wrote 90%'
               do y=1,nyp
                  do x=1,nx/2
                     urx(2*x-1)=boxr(x,z-zb+1,y)
                     urx(2*x)=boxi(x,z-zb+1,y)
                  end do
                  write(11) urx
               end do
            end do
         end do
      end do
      close(unit=11)

      return

      end subroutine wdisca
