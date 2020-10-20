c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine blvel(bu,xsc,x0,eta,dybla,dfbla,gbl,th,nbla,rlam,
     &     spanv,re,xl,bstart,blength,u0low,prex,prey,mx,my,wb,
     &     xbl1,xbl2,bst,bc,m1,scalar,thgrad,pr)

      implicit none

      integer nbla,mx,my,scalar
      real th(nbla,scalar),dfbla(nbla),gbl(nbla),dybla
      real x0,re,xl,xsc,bstart,blength,u0low,rlam,spanv
      real bu(mx+2,my,3+scalar)
      real eta(my)
      real prex(mx+15),prey(my*2+15)
      real wb(mx+2,my)

      integer x,y,j
      real xbl1(mx),xbl2(mx)
      real xbl,xmin,uu1,uu2,thuu1,thuu2
      real bst(mx)
      real xc,ybl,etabl,df1,df2,tmp
      real bc(mx)
      real m1(scalar),thgrad(scalar),thgrad1,pr(scalar)
      real t1,t2
c
c     Functions
c
      real cubip,step,xlim
c
c     Constant
c
      real pi
      parameter (pi = 3.1415926535897932385)
c
c     Calculate base flow in bu(x,y,1..3)
c     Calculate the velocity profile for u
c     by blending the flow from the upstream
c     and the downstream part we need a coordinate which starts at
c     bstart and goes on to fend+xl
c
c     To avoid complex looping we first construct the x-coordinates
c     for each point since we need to blend the 'upstream' and
c     'downstream' boundary layer we construct one coordinate which
c     is one period higher
c
c     bst is the blending step function
c     xbl1 goes from x0+bstart to x0+bstart+xl
c     xbl2 goes from x0+bstart+xl to x0+bstart+2xl
c
      do x=1,mx
         xc=real(x-mx/2-1)/real(mx)*xl+xsc
         xbl=xc-(int((xc-bstart)/xl+1.)-1.)*xl
         bst(x)=step((xbl-bstart)/blength)
         xbl1(x)=xbl+x0
         xbl2(x)=xbl1(x)+xl
      end do
c
c     Now limit the coordinate to avoid going upstream of leading edge
c
      xmin=x0/2.
      do x=1,mx
         xbl1(x)=xlim(xbl1(x),xmin)
         xbl2(x)=xlim(xbl2(x),xmin)
      end do
c
c     Now compute u,w,th from the similarity solution
c
      do y=1,my
         ybl=1.+eta(y)
         do x=1,mx
c
c     Schlichting scaled eta from xbl1
c
            uu1=(xbl1(x)/x0)**rlam
            etabl=ybl/xbl1(x)*sqrt((rlam+1.)*xbl1(x)*uu1*re/2.)
c
c     Interpolate u
c
            df1=cubip(etabl,dfbla,dybla,nbla)
c
c     Schlichting scaled eta from xbl2
c
            uu2=(xbl2(x)/x0)**rlam
            etabl=ybl/xbl2(x)*sqrt((rlam+1.)*xbl2(x)*uu2*re/2.)
c
c     Interpolate u
c
            df2=cubip(etabl,dfbla,dybla,nbla)
c
c     Compute base flow
c
            bu(x,y,1)=uu1*df1*bst(x)+uu2*df2*(1.-bst(x))+u0low
c
c     Schlichting scaled eta from xbl1
c
            etabl=ybl/xbl1(x)*sqrt((rlam+1.)*xbl1(x)*uu1*re/2.)
c
c     Interpolate w
c
            df1=cubip(etabl,gbl,dybla,nbla)
c
c     Schlichting scaled eta from xbl1
c
            etabl=ybl/xbl2(x)*sqrt((rlam+1.)*xbl2(x)*uu2*re/2.)
c
c     Interpolate w
c
            df2=cubip(etabl,gbl,dybla,nbla)
            bu(x,y,3)=spanv*(df1*bst(x)+df2*(1.-bst(x)))

            do j=1,scalar
c
c     Schlichting scaled eta from xbl1
c
               thuu1=(xbl1(x)/x0)**m1(j)
               etabl=ybl/xbl1(x)*sqrt((rlam+1.)*xbl1(x)*uu1*re/2.)
c
c     Interpolate th
c
               df1=cubip(etabl,th(1,j),dybla,nbla)
c
c     Schlichting scaled eta from xbl1
c
               thuu2=(xbl2(x)/x0)**m1(j)
               etabl=ybl/xbl2(x)*sqrt((rlam+1.)*xbl2(x)*uu2*re/2.)
c
c     Interpolate w
c
               df2=cubip(etabl,th(1,j),dybla,nbla)
               bu(x,y,3+j)=thuu1*df1*bst(x)+thuu2*df2*(1.-bst(x))
            end do

         end do
      end do
c
c     Now construct the v velocity from u by continuity in Chebyshev space
c
      do y=1,my
         do x=1,mx
            bu(x,y,2) = bu(x,y,1)
         end do
      end do
c
c     First fourier transform
c
      call vrfftf(bu(1,1,2),bu(2,1,2),wb,wb(2,1),mx,my,2,mx+2,prex)
c
c     Set odd-ball to zero
c
      do y=1,my
         bu(mx+1,y,2)=0.0
      end do
c
c     Chebyshev transform
c
      call vchbf(bu(1,1,2),wb,my,mx,mx+2,1,prey)
c
c     Normalize
c
      uu1=2./(real(mx)*real(my-1))
      do y=1,my
         do x=1,mx
            bu(x,y,2)=bu(x,y,2)*uu1
         end do
      end do
c
c     Compute dvdy = -dudx = -i alpha u
c     It's a simple first-order ODE with known RHS.
c
      do y=1,my
         do x=1,mx/2
            tmp         =2.*pi/xl*real(x-1)*bu(2*x,y,2)
            bu(2*x,y,2)=-2.*pi/xl*real(x-1)*bu(2*x-1,y,2)
            bu(2*x-1,y,2)=tmp
         end do
      end do
c
c     Then integrate
c     bc(x) is the physical boundary condition, v=0 at y=0.
c
      do x=1,mx
         bc(x)=0.0
      end do
      call icheb(wb,bu(1,1,2),bc,my,mx,mx+2)
      do y=1,my
         do x=1,mx
            bu(x,y,2)=wb(x,y)
         end do
      end do
c
c     Now find the value at the lower boundary for each mode
c     note skipping of mode 1
c
      do y=2,my-1,2
         do x=1,mx
            bc(x)=bc(x)+bu(x,y,2)-bu(x,y+1,2)
         end do
      end do
c
c     To set bc v=0 for y=0
c
      do x=1,mx
         bu(x,1,2)=bc(x)
      end do
c
c     Go back to physical space
c     Chebyshev transform
c
      call vchbb(bu(1,1,2),wb,my,mx,mx+2,1,prey)
c
c     Fourier transform
c
      call vrfftb(bu(1,1,2),bu(2,1,2),wb,wb(2,1),mx,my,2,mx+2,prex)
c
c     Rescaling of scalar field for constant wall flux
c
      do j=1,scalar
         write(*,*) 'Scalar ',j,' Pr=',pr(j),' m1=',m1(j)
         if (abs(m1(j)-0.5).lt.1e-13) then
            write(*,*) 'm1=0.5 detected --> constant wall gradient'
            read(10,*) thgrad1
            read(10,*) t1
            write(*,*) 'Rescale theta to prescribed wall gradient'
            write(*,*) '  thgrad (old) = ',thgrad(j)
            write(*,*) '  thgrad (new) = ',thgrad1
            write(*,*) 'freestream val.= ',t1
            do y=1,my
               do x=1,mx
                  bu(x,y,3+j) = t1+bu(x,y,3+j) * (thgrad1/thgrad(j))
               end do
            end do
         end if
         if (m1(j).eq.0) then
            write(*,*) 'm1=0 detected'
            read(10,*) t1
            read(10,*) t2
            write(*,*) 'wall value       : ',t1
            write(*,*) 'freestream value : ',t2
            do y=1,my
               do x=1,mx
                  bu(x,y,3+j) = t2+(t1-t2)*bu(x,y,3+j)
               end do
            end do
         end if
      end do

      end subroutine blvel
