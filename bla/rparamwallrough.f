c ***********************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      subroutine rparamwallrough(wallroughinput,updat,every,
     &     v_wall,taylor4,h_rough,hstart,hend,hrise,hfall,hrms,
     &     zrand,rseed,nzr,psi,rghfil,roughfile,pert,my_node)
c
c     Reads in the parameters needed for the modification
c     (projection based on Taylor expansion) of the boundary
c     conditions at the "wall", if there is surface roughness.
c
c     The file which to read from is specified in rparambl.f.
c
c     Called by bla.f
c
      implicit none

      include 'par.f'

      integer every,nzr,rseed,my_node
      real h_rough,hstart,hend,hrise,hfall,psi
      logical updat,v_wall,taylor4,rghfil,hrms,zrand,pert
      character(32) wallroughinput,roughfile
c
c     Initialization
c
      updat=.false.
      every=0
      v_wall=.false.
      taylor4=.false.
      h_rough=0.
      hstart=0.
      hend=0.
      hrise=0.
      hfall=0.
      hrms=.false.
      zrand=.false.
      nzr = 1
      rseed=-1
      psi=0.
      rghfil=.false.


      if (my_node.eq.0) then
         write(ios,*)
         write(ios,*) '--> Wall roughness included. Opening  ',
     &        trim(wallroughinput)
         write(ios,*)
      end if
c
c     Read the wall-roughness parameters
c
      open(unit=15,file=wallroughinput,status='old')

      read(15,*) updat

      if (updat) then
         if (pert) then
            updat=.false.
            if (my_node .eq. 0) then
               write(ios,*)
               write(ios,*) '--> Note: Update of roughness conditions'
               write(ios,*) '    not possible in perturbation mode'
               write(ios,*)
            end if
         end if

         read(15,*) every
         if (.not.pert) then
            if (my_node.eq.0) then
               write(ios,*)
               write(ios,*) 'Wall-roughn. BCs updated every ',
     &              every,'th dt'
               write(ios,*)
            end if
         end if
      end if

      read(15,*) v_wall
      read(15,*) taylor4
      read(15,*) h_rough
      read(15,*) hstart
      read(15,*) hend
      read(15,*) hrise
      read(15,*) hfall
      read(15,*) hrms
      read(15,*) zrand

      if (zrand) read(15,*) rseed

      read(15,*) nzr
      read(15,*) psi

      read(15,*) rghfil
      if (rghfil) read(15,*) roughfile

      close(15)


      end subroutine  rparamwallrough
