c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine statpl(plxy,plmxy,plrxy,plxz,plmxz,plrxz,
     &plyz,plmyz,plryz,nzn,tint,nint,tf,tl,ipln,tpl,sym)
c
c     Calculate statistics from plane sequences
c     this routine recieves one plane at the time and accumulates
c     average and squared average
c     at the last plane rms is calculated
c
      implicit none

      include 'par.f'

      integer nzn,nint,ipln,tpl
      real tint(nint),tf,tl,sym
      real plxy(nx+1,nyp),plmxy(nx+1,nyp),plrxy(nx+1,nyp)
      real plxz(nx+1,nz+1),plmxz(nx+1,nz+1),plrxz(nx+1,nz+1)
      real plyz(nyp,nz+1),plmyz(nyp,nz+1),plryz(nyp,nz+1)

      real c
      integer x,y,z
c
c     Calculate integration weights
c
      c=1.
      call intwgt(c,tint,nint,tf,tl,ipln)
c
c     Reset statistics arrays
c
      if (ipln.eq.1) then
         if (tpl.eq.1) then
            do y=1,nyp
               do x=1,nx+1
                  plmxy(x,y)=0.
                  plrxy(x,y)=0.
               end do
            end do
         elseif (tpl.eq.2) then
            do z=1,nzn
               do x=1,nx+1
                  plmxz(x,z)=0.
                  plrxz(x,z)=0.
               end do
            end do
         else
            do y=1,nyp
               do z=1,nzn
                  plmyz(y,z)=0.
                  plryz(y,z)=0.
               end do
            end do
         end if
      end if
c
c     Add upp statistics
c
c      write(*,*) c,nint,ipln
      if (tpl.eq.1) then

         do y=1,nyp
            do x=1,nx+1
               plmxy(x,y)=plmxy(x,y)+c*plxy(x,y)
               plrxy(x,y)=plrxy(x,y)+c*plxy(x,y)**2
            end do
         end do
      elseif (tpl.eq.2) then
         do z=1,nzn
            do x=1,nx+1
               plmxz(x,z)=plmxz(x,z)+c*plxz(x,z)
               plrxz(x,z)=plrxz(x,z)+c*plxz(x,z)**2
            end do
         end do
      else
         do y=1,nyp
            do z=1,nzn
               plmyz(y,z)=plmyz(y,z)+c*plyz(y,z)
               plryz(y,z)=plryz(y,z)+c*plyz(y,z)**2
            end do
         end do
      end if
c
c     Last plane calculate rms
c
      if (ipln.eq.nint) then
         if (tpl.eq.1) then
            do y=1,nyp
               do x=1,nx+1
cc            if (plrxy(x,y)-plmxy(x,y)**2.lt.0)
cc     &   write(*,*) x,y,ipln,plrxy(x,y),plmxy(x,y),
cc     &plrxy(x,y)-plmxy(x,y)**2
                  plrxy(x,y)=sqrt(max(0.,plrxy(x,y)-plmxy(x,y)**2))
               end do
            end do
         elseif (tpl.eq.2) then
            do z=1,nzn
               do x=1,nx+1
                  plrxz(x,z)=sqrt(max(0.,plrxz(x,z)-plmxz(x,z)**2))
               end do
            end do
            if (nfzsym.eq.1) then
               do z=nz/2+2,nz
                  do x=1,nx+1
                     plmxz(x,z)=plmxz(x,nz+2-z)*sym
                     plrxz(x,z)=plrxz(x,nz+2-z)*sym
                  end do
               end do
            end if
c          do z=1,nzn
c             do x=1,nx+1
c                plmxz(x,z) = plmxz(x,z)/nint
c             end do
c          end do
            do x=1,nx+1
               plmxz(x,nz+1)=plmxz(x,1)
               plrxz(x,nz+1)=plrxz(x,1)
            end do
         else
            do y=1,nyp
               do z=1,nzn
                  plryz(y,z)=sqrt(max(0.,plryz(y,z)-plmyz(y,z)**2))
               end do
            end do
            if (nfzsym.eq.1) then
               do y=1,nyp
                  do z=nz/2+2,nz
                     plmyz(y,z)=plmyz(y,nz+2-z)*sym
                     plryz(y,z)=plryz(y,nz+2-z)*sym
                  end do
               end do
            end if
            do y=1,nyp
               plmyz(y,nz+1)=plmyz(y,1)
               plryz(y,nz+1)=plryz(y,1)
            end do
         end if
      end if

      return

      end subroutine statpl
