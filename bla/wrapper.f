C************************************************************************
c
c $HeadURL$
c $LastChangedDate$
c $LastChangedBy$
c $LastChangedRevision$
c
c ***********************************************************************
      program wrapper


      implicit none
      
#include "par.f"
#ifdef MPI
      include 'mpif.h'
#endif

c
c     Main storage (distributed in z among the processors)
c
      real ur(memnx,memny,memnz,memnxyz),ui(memnx,memny,memnz,memnxyz)
c
c     Variables for preprbl
c
      real alfa(nx/2,mbz),beta(nz)
      real eta(nyp),deta(nyp),wint(nyp),wd1(nyp,4),xl,zl
      real gridx(nx),gridy(nyp),gridz(nz)
      real prex(nxp+15),prey(nyp*2+15)
      real prez(nzp*2+15),pres(nzst*2+15),prea(nzat*3/2+15)
      real prexn(nx+15),prezn(nz*2+15),presn(nz+2+15),prean(nz*3/4+15)
      real prezr(nzp+15)
      real wr  ((nxp/2+1)*mby*nzd         ,nthread)
      real wi  ((nxp/2+1)*mby*nzd         ,nthread)
      logical gall
      real deltaxyz2(nyp)
c
c     Variables for rdiscbl
c
      character*80 w_namin,w_namout,b_namin,startname
      real re,t
      real u0low,u0upp,w0low,w0upp,du0upp
      real zs,xs,dstar,bstart,bslope,spanv,rlam,pr(scalar),m1(scalar)
      logical spat,varsiz,pert
      real boxr(nx/2,mbz,nyp),boxi(nx/2,mbz,nyp)
      real w3(nx/2*mbz*nyp,nthread),urx(nx)
      integer nxtmp,fltype
      real br(memnx,memny,3),bi(memnx,memny,3)
      logical rbflow

c
c     Wrapper variables
c
      real tmin,alam,aturb
      real acc,acc_lim,lambda,lambdaLM,lambdaTM,lambdaTM0
      logical tstatus
c
c     Extra variables
c
      integer x,y,z,i
      character*20 sttimec

c     ################################################################
c
c     MPI
c
      integer all_nodes,my_node,ierror

c     ################################################################

c-----------------------------------------------------------------------
c     BEGINNING OF PROGRAM
c-----------------------------------------------------------------------
c
c     Initialize timer
c
      call time_string(sttimec)

#ifdef MPI
c
c     Startup MPI
c
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,my_node,ierror)
      call mpi_comm_size(mpi_comm_world,all_nodes,ierror)
#else
c
c     Set MPI variables to serial settings
c
      my_node = 0
      all_nodes = 1
#endif

      if (my_node.eq.0) then
         write(ios,*) '********************************************'//
     &              '***********************'
         write(ios,*) '*                                           '//
     &              '                      *'
         write(ios,*) '*                             Simson        '//
     &              '                      *'
         write(ios,*) '*               wrapper for bla $Rev$'//
     &              '                      *'
         write(ios,*) '*                                           '//
     &              '                      *'
         write(ios,*) '********************************************'//
     &              '***********************'

         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  General information  <<<<<<<'
         write(ios,*) '--------------------------------------------'//
     &        '-----------------------'
         write(ios,*) 'Started              : ',sttimec
      end if

c
c     Reading wrapper.i file
c
      open(unit=101,status='old',file='wrapper.i')

      if (my_node.eq.0) then
         write(ios,*)
         write(ios,*)
         write(ios,*) '>>>>>>  Parameters read from wrapper.i <<<<<<<'
         write(ios,*) '--------------------------------------------'//
     &              '-----------------------'
      end if

      call comment(101)
      read(101,*) w_namin
      if (my_node.eq.0) then
         write(ios,*) 'Initial field (w_namin)         : ',trim(w_namin)
      end if
      
      call comment(101)
      read(101,*) b_namin    
      if (my_node.eq.0) then
         write(ios,*) 'Base field (b_namin)            : ',trim(b_namin)
      end if
      
      call comment(101)
      read(101,*) w_namout 
      if (my_node.eq.0) then
         write(ios,*) 'End field (w_namout)            : ',
     &        trim(w_namout)
      end if
      
      call comment(101)
      read(101,*) spat
      if (my_node.eq.0) then
         write(ios,*) 'Spatial simulation (spat)       :',spat
      end if

      call comment(101)
      read(101,*) pert
      if (my_node.eq.0) then
         write(ios,*) 'Perturbation mode (pert)        :',pert
      end if

c      
c     Bissection loop variables (e.g. for restarting purposes)
c
      call comment(101)
      read(101,*) tmin
      if (my_node.eq.0) then
         write(ios,*) 'Minimum check time (tmin)       :',tmin
      end if

      call comment(101)
      read(101,*) alam
      if (my_node.eq.0) then
         write(ios,*) 'Laminar amplitude (alam)        :',alam
      end if

      call comment(101)
      read(101,*) aturb
      if (my_node.eq.0) then
         write(ios,*) 'Turbulent amplitude (aturb)     :',aturb
      end if

      call comment(101)
      read(101,*) lambda
      if (my_node.eq.0) then
         write(ios,*) 'Initial multiplier (lambda)     :',lambda
      end if
     
      call comment(101) 
      read(101,*) lambdaLM
      if (my_node.eq.0) then
         write(ios,*) 'Laminar multiplier (lambdaLM)   :',lambdaLM
      end if
      
      call comment(101)
      read(101,*) lambdaTM0
      if (my_node.eq.0) then
         write(ios,*) 'Turbulent multiplier (lambdaTM) :',lambdaTM0
      end if

      call comment(101)
      read(101,*) acc_lim
      if (my_node.eq.0) then
         write(ios,*) 'Accuracy (acc_lim)              :',acc_lim
      end if

      call comment(101)
      read(101,*) rbflow
      if (my_node.eq.0) then
         write(ios,*) 'Rescale using base flow (rbflow):',rbflow
      end if
      
      close(unit=101)
            
c     ################################################################

      t = 0.
      acc = abs(lambdaTM-lambdaLM)
      lambdaTM = lambdaTM0

      if (rbflow) then
c
c     Read base flow and store relevant plane
c     We consider a baseflow which is a function of x and y.
c
         call rdiscbl(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,pr,
     &        xl,zl,t,xs,dstar,fltype,bstart,bslope,rlam,m1,spanv,
     &        spat,b_namin,varsiz,3+scalar,boxr,boxi,w3,urx,
     &        nxtmp,my_node)      
         if (my_node.eq.0) then
            do i=1,3
               br(:,:,i) = ur(:,:,1,i)
               bi(:,:,i) = ui(:,:,1,i)
            end do
         else
            br = 0.
            bi = 0.
         end if
      else
         br = 0.
         bi = 0.
      end if

c
c     Write initial information
c
      if (my_node.eq.0) then
         write(ios,*)  'lambda',
     &        lambda,lambdaLM,lambdaTM,t/dstar
         open(file='lambda.out',unit=103)
         write(103,'(4e25.16,f18.9)')  lambda,lambdaLM,lambdaTM,acc,
     &        t/dstar
         close(103)
      end if

c
c     Main wrapper loop
c
      do while(acc.gt.acc_lim)
c
c     Read initial velocity field
c
         call rdiscbl(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,pr,
     &        xl,zl,t,xs,dstar,fltype,bstart,bslope,rlam,m1,spanv,
     &        spat,w_namin,varsiz,3+scalar,boxr,boxi,w3,urx,
     &        nxtmp,my_node)      

c
c     Reset time to zero
c
         t = 0.
c
c     Initialize certain preprocessing information for FFTs, grid etc.
c
         call preprbl(alfa,beta,eta,deta,xl,zl,prex,prey,prez,pres,prea,
     &        prexn,prezn,presn,prean,prezr,wint,wd1,
     &        gridx,gridy,gridz,dstar,deltaxyz2)
c
c     Scale disturbance
c
         if (my_node.eq.0) then
c     
c     Processor which includes zero/zero mode
c
            if (rbflow) then
c
c     Subtract baseflow, rescale, add baseflow, ....
c     This is done for all wavenumbers with beta=0.
c
               do i=1,3
                  do z=1,1
                     do y=1,memny
                        do x=1,memnx
                           ur(x,y,z,i) = br(x,y,i) +
     &                          lambda*(ur(x,y,z,i) - br(x,y,i))
                           ui(x,y,z,i) = bi(x,y,i) +
     &                          lambda*(ui(x,y,z,i) - bi(x,y,i))
                        end do
                     end do
                  end do
c
c     Only rescale for beta.ne.0
c
                  do z=2,memnz
                     do y=1,memny
                        do x=1,memnx
                           ur(x,y,z,i)=lambda*ur(x,y,z,i)
                           ui(x,y,z,i)=lambda*ui(x,y,z,i)
                        end do
                     end do
                  end do
               end do
            else
c
c     ..., or rescale all except 0/0 mode
c
               
               do i=1,3
                  do z=1,1
                     do y=1,memny
                        do x=2,memnx
                           ur(x,y,z,i)=lambda*ur(x,y,z,i)
                           ui(x,y,z,i)=lambda*ui(x,y,z,i)
                        end do
                     end do
                  end do
                  do z=2,memnz
                     do y=1,memny
                        do x=1,memnx
                           ur(x,y,z,i)=lambda*ur(x,y,z,i)
                           ui(x,y,z,i)=lambda*ui(x,y,z,i)
                        end do
                     end do
                  end do
               end do
            end if
         else
c
c     Processor without zero/zero mode
c     In the 1D parallelisation, the mean flow is always on my_node=0.
c
            do i=1,3
               do z=1,memnz
                  do y=1,memny
                     do x=1,memnx
                        ur(x,y,z,i)=lambda*ur(x,y,z,i)
                        ui(x,y,z,i)=lambda*ui(x,y,z,i)
                     end do
                  end do
               end do
            end do
         end if
            


c
c     Write out scaled infile
c         
         startname = "in.u"
         call wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &        bstart,bslope,rlam,spanv,startname,3+scalar,gall,
     &        boxr,boxi,urx,alfa,zs,beta,my_node)

         call bla(ur,ui,u0low,u0upp,w0low,w0upp,du0upp,re,pr,
     &        xl,zl,t,xs,dstar,fltype,bstart,bslope,rlam,m1,spanv,
     &        spat,varsiz,w3,urx,nxtmp,my_node,all_nodes,
     &        tstatus,tmin,alam,aturb)
         
         call wdiscbl(ur,ui,re,pr,m1,xl,zl,t,xs,dstar,fltype,
     &        bstart,bslope,rlam,spanv,w_namout,3+scalar,gall,
     &        boxr,boxi,urx,alfa,zs,beta,my_node)
        
         if (my_node.eq.0) then
c
c       loop for bissection--------------------------------------------
c
            if (tstatus) lambdaTM = min(lambda,lambdaTM)
            if (.not.tstatus) lambdaLM = max(lambda,lambdaLM)
            acc = abs(lambdaTM-lambdaLM)
            
            if (lambdaTM.le.lambdaTM0) then
               lambda = (lambdaTM+lambdaLM)*0.5
            else
               lambda = 2*lambdaLM
            end if
         end if
c         
c     Communicate lambda to all nodes
c        
#ifdef MPI
         if (nproc.gt.1) then
            call mpi_bcast(lambda,1,mpi_double_precision,0,
     &           mpi_comm_world,ierror)
            call mpi_bcast(acc,1,mpi_double_precision,0,mpi_comm_world,
     &           ierror)
         end if
#endif

         if (my_node.eq.0) then
            open(file='wrapper.out',unit=102,status='replace')
            write(102,*) w_namin
            write(102,*) w_namout
            write(102,*) spat
            write(102,*) pert
            write(102,*) tmin
            write(102,*) alam
            write(102,*) aturb
            write(102,*) lambda
            write(102,*) lambdaLM
            write(102,*) lambdaTM
            write(102,*) acc_lim
            write(102,*) rbflow
            close(102)
         end if
         
c
c     Write current step information
c
         if (my_node.eq.0) then
            write(ios,*)  'lambda',
     &           lambda,lambdaLM,lambdaTM,t/dstar
            open(file='lambda.out',unit=103,access='append')
            write(103,'(4e25.16,f18.9)')  lambda,lambdaLM,lambdaTM,acc,
     &           t/dstar
            close(103)
         end if  
      end do

c     ################################################################

c
c     Finalize program
c
#ifdef MPI
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_finalize(ierror)
#endif
      end program wrapper
