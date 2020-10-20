c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine intwgt(c,tint,nint,tf,tl,ipln)
c
c     Calculate integration weight for the plane ipln
c
      implicit none

      real c
      integer nint,ipln
      real tint(nint),tf,tl

      if (nint.eq.1) c=1.
      if (nint.eq.2) then
         if (ipln.eq.1) c=(2.*tint(2)-tf-tl)*.5/(tint(2)-tint(1))
         if (ipln.eq.2) c=-(2.*tint(1)-tf-tl)*.5/(tint(2)-tint(1))
      end if
      if (nint.eq.3) then
         if (ipln.eq.1) c=(tint(2)-tf)**2/(tint(2)-tint(1))*.5/(tl-tf)
         if (ipln.eq.2) c=
     &        ((tint(2)-tf)*(tint(2)-2.*tint(1)+tf)/(tint(2)-tint(1))+
     &        (tint(2)-tl)*(tint(2)-2.*tint(3)+tl)/(tint(3)-tint(2)))*
     &        .5/(tl-tf)
         if (ipln.eq.3) c=(tint(2)-tl)**2/(tint(3)-tint(2))*.5/(tl-tf)
      end if

      if (nint.gt.3) then

c      if (ipln.eq.1) c=(tint(2)-tf)**2/(tint(2)-tint(1))*.5/(tl-tf)

c      if (ipln.eq.2) c=((tint(2)-tf)*(tint(2)-2.*tint(1)+tf)/
c     &(tint(2)-tint(1))+tint(3)-tint(2))*.5/(tl-tf)

c      if (ipln.gt.2.and.ipln.lt.nint-1)
c     &     c=(tint(ipln+1)-tint(ipln-1))*.5/(tl-tf)


c      if (ipln.eq.nint-1)
c     &c=-((tint(nint-1)-tl)*(tint(nint-1)-2.*tint(nint)+tl)/
c     &(tint(nint-1)-tint(nint))+tint(nint-2)-tint(nint-1))*.5/(tl-tf)

c      if (ipln.eq.nint)
c     &c=-(tint(nint-1)-tl)**2/(tint(nint-1)-tint(nint))*.5/(tl-tf)


         if (ipln.gt.1.and.ipln.lt.nint)
     &        c= .5*(tint(ipln+1)-tint(ipln-1))/(tint(nint)-tint(1))
         if (ipln.eq.1)
     &        c=    (tint(ipln+1)-tint(ipln))  /(tint(nint)-tint(1))
         if (ipln.eq.nint)
     &        c=    (tint(ipln)-tint(ipln-1))  /(tint(nint)-tint(1))

      end if

      return

      end subroutine intwgt
