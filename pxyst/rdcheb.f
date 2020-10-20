c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
c
      subroutine rdcheb(a,n,m,md)
c      calculates first derivative in chebychev space
c
c    belongs to pxxys version 1.0
c    for more info see the pxxys.f file
c
c      n number of chebychev coefficients
c      m number of derivatives
c      md first index length on a array
c
      implicit logical (a-z)
      integer n,m,md
      real a(md,n)
c
      integer i,k
      real tmp1
      do 2000 k=1,m
        tmp1=a(k,n-2)
        a(k,n-2)=2.*real(n-2)*a(k,n-1)
        a(k,n-1)=2.*real(n-1)*a(k,n)
        a(k,n)=tmp1
 2000 continue
      do 2010 i=n-3,1,-1
        do 2020 k=1,m
          tmp1=a(k,n)
          a(k,n)=a(k,i)
          a(k,i)=a(k,i+2)+2.*real(i)*tmp1
 2020   continue
 2010 continue
      do 2030 k=1,m
         a(k,n)=0.0
         a(k,1)=.5*a(k,1)
 2030 continue
      return
      end
