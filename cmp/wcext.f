c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine wcext(vext,cext,u0low,dstar,eta,fltype,xl)
c
c     Print out the maximum of velocities and vorticities and their coordinates
c
      implicit none

      include 'par.f'

      integer fltype
      real vext(nyp,2,6),cext(nyp,2,6,2)
      real dstar,eta(nyp),u0low,xl
c
      integer y,i,j,k
      real vextf(2,6),cextf(2,6,3),vm
      character*32 form1
      character*80 form

      do i=1,6
         vextf(1,i)=1.E20
         vextf(2,i)=-1.E20
      end do
      do y=1,nyp
         do i=1,6
            vm=0.0
            if (i.eq.1) then
               if (fltype.eq.1.or.fltype.eq.4) vm=u0low+1.-eta(y)**2
               if (fltype.eq.2.or.fltype.eq.5) vm=eta(y)
c
c     Mean values not implemented for bl-flows
c
               if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) vm=u0low
            end if
            if (i.eq.6) then
               if (fltype.eq.1.or.fltype.eq.4) vm=2.*eta(y)
               if (fltype.eq.2.or.fltype.eq.5) vm=-1.
c
c     Mean values not implemented for bl-flows
c
               if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) vm=0.0
            end if
            if (vextf(1,i).gt.vext(y,1,i)-vm) then
               vextf(1,i)=vext(y,1,i)-vm
               cextf(1,i,1)=cext(y,1,i,1)
               cextf(1,i,2)=eta(y)
               cextf(1,i,3)=cext(y,1,i,2)
            end if
            if (vextf(2,i).lt.vext(y,2,i)-vm) then
               vextf(2,i)=vext(y,2,i)-vm
               cextf(2,i,1)=cext(y,2,i,1)
               cextf(2,i,2)=eta(y)
               cextf(2,i,3)=cext(y,2,i,2)
            end if
         end do
      end do
c
c     Spatial simulations the x-coordinate in (0,xl)
c
      if (fltype.ge.4) then
         do i=1,2
            do j=1,6
               cextf(i,j,1)=mod(cextf(i,j,1)+xl,xl)
            end do
         end do
      end if
c
      write(*,*) 'Extremum amplitudes'
      do i=1,2
         do j=1,6
c
c     omx and omz not implemented yet
c
            if (i.eq.1.and.j.eq.1) form1='(''min of u at x,y,z  '''
            if (i.eq.1.and.j.eq.2) form1='(''min of v at x,y,z  '''
            if (i.eq.1.and.j.eq.3) form1='(''min of w at x,y,z  '''
            if (i.eq.1.and.j.eq.4) form1='(''min of omx at x,y,z'''
            if (i.eq.1.and.j.eq.5) form1='(''min of omy at x,y,z'''
            if (i.eq.1.and.j.eq.6) form1='(''min of omz at x,y,z'''
            if (i.eq.2.and.j.eq.1) form1='(''max of u at x,y,z  '''
            if (i.eq.2.and.j.eq.2) form1='(''max of v at x,y,z  '''
            if (i.eq.2.and.j.eq.3) form1='(''max of w at x,y,z  '''
            if (i.eq.2.and.j.eq.4) form1='(''max of omx at x,y,z'''
            if (i.eq.2.and.j.eq.5) form1='(''max of omy at x,y,z'''
            if (i.eq.2.and.j.eq.6) form1='(''max of omz at x,y,z'''
            form=form1//',F20.14,3F10.4)'
            if (fltype.eq.3.or.fltype.ge.6.or.fltype.lt.0) then
               if (j.le.3) then
                  write(*,form) vextf(i,j),cextf(i,j,1)/dstar,
     &                 (1.+cextf(i,j,2))/dstar,cextf(i,j,3)/dstar
               else
                  write(*,form) vextf(i,j)*dstar,cextf(i,j,1)/dstar,
     &                 (1.+cextf(i,j,2))/dstar,cextf(i,j,3)/dstar
               end if
            else
               write(*,form) vextf(i,j),(cextf(i,j,k),k=1,3)
            end if
         end do
      end do

      return

      end subroutine wcext
