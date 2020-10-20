c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      real function ran2(idum)
c
c     A simple portable random number generator
c
c     Requires 32-bit integer arithmetic
c     Taken from Numerical Recipes, William Press et al.
c     gives correlation free random numbers but does not have a very large
c     dynamic range, i.e only generates 714025 different numbers
c     for other use consult the above
c     Set idum negative for initialization
c
      implicit none

      integer idum,j

      integer,parameter :: m=714025,ia=1366,ic=150889
      real,parameter :: rm=1./m

      integer,save :: ir(97),iy
      integer,save :: iff = 0

      if (idum.lt.0.or.iff.eq.0) then
c
c     Initialize
c
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
c
c     Generate random number
c
      j=1+(97*iy)/m
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum

      end function ran2
