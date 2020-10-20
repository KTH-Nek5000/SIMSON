#!/bin/tcsh
# The bla.i needs to be as follows:
# initial field: start.u or start.uu
# end field: end.u or end.uu
# statistics: xy.stat
# correlations: corr.stat and corrx.stat
# time series: ser.stat
# amplitudes: amp.stat
#
# Philipp Schlatter, 2007
#
#------------------------------------------------

# move to the appropriate directory
cd /cfs/ekman/scratch/p/pschlatt/geert/run_768x385x768

# indicates that the script is running
echo RUNNING > STATUS

# output is appended to script.out
echo BLA RUNSCRIPT $SP_JID
echo -------------------------
pwd  
date 
echo -------------------------

set i = 1000
while (${i} >= 1)
  if ( -e inter.uu.${i}-00000 ) then
    break
  endif  
  set i = "`expr ${i} - 1`"
end
echo THIS IS ${i}


# Presently, 1000 is the maximum number of fields. This is set to
# avoid redoing over and over calculation 1000.
if ( $i > 999) then
	echo Nothing to do. Job $i is too big.
	echo OVERFLOW > STATUS
        exit 1
else

# Check if the initial file exists. Otherwise run bls or copy init.u
        if ( -e inter.uu.${i}-01024 ) then
        else
          echo ERROR > STATUS
          exit 1
          if ( -e init.uu ) then
	    echo using init.uu
          else
	    echo First do bls
            ./bls
          endif
          ./blamv.sh init.uu inter.uu.0
        endif

# Prepare the run by deleting old working files
	rm -f start.uu-?????
	./blaln.sh inter.uu.${i} start.uu
	rm -f history.out
	rm -f step.out
	rm -f end.uu-?????
	rm -f end.uu.p

	echo Starting job from inter.uu.$i
	echo RUNNING from $i > STATUS

# Run the job. Here use the appropriate command for the architecture
# this is for Lucidor:
#      ./mpich.lxl -p 2 -l ./bla.lucidor	
# this is for Lenngren:
#     /opt/scali/bin/mpirun -np 64 -npn 2 -inherit_limits ./bla.lenngren > bla.out
# this is for serial
#        ./bla > bla.out
# this is for the BG

    module add beta-modules
    module add pgi/10.3
    module add openmpi/1.4.2-pgi

    cp -f $SP_HOSTFILE old_host_file.$SP_JID

    bash -c 'cat old_host_file.$SP_JID | while read line; do for ((i=1;i<=8;i+=1)); do echo "${line}"; done; done > new_host_file.$SP_JID'

    mpirun --mca mpi_leave_pinned 0 --mca mpi_paffinity_alone 1 -n 1024  -machinefile ./new_host_file.$SP_JID ./blaq.1024 > bla.out.$SP_JID





# Job has finished
	echo -- job $i finished

	set next = "`expr ${i} + 1`"
	set prev = "`expr ${i} - 1`"

# Copy output files
	mv -f bla.out.$SP_JID bla.out.${next}
	if( -e end.uu-01024 ) then
    		echo -- output from job $i detected

                ./blamv.sh  end.uu inter.uu.${next}
		if( -e end.uu.p ) then
		       mv -f end.uu.p inter.uu.p.${next}
                endif

		if( -e inter.uu.${prev}-00000 ) then
                       rm -f inter.uu.${prev}-?????
                endif

		if( -e xy.stat ) then
	                mv -f xy.stat xy.stat.${next}
		endif
		if( -e corr.stat ) then
			mv -f corr.stat corr.stat.${next}
		endif
		if( -e corrx.stat ) then
			mv -f corrx.stat corrx.stat.${next}
		endif
		if( -e ser.stat) then
			mv -f ser.stat ser.stat.${next}
		endif
		if( -e amp.stat) then
			mv -f amp.stat amp.stat.${next}
		endif
#                if( -e spec.stat-00000) then
#                        mkdir spec.${next}
#                        mv -f spec.stat-????? spec.${next}
#                        mv -f spec_.stat-????? spec.${next}
#                endif

                mv -f step.out step.out.${next}
                mv -f history.out history.out.${next}
                rm -f start.uu-?????

   	        echo -- output copied in XXX.${next}

# Here you could add a resubmit command (not needed on Lenngren)
                if( -e stop.now ) then
		   echo STOPPED job $i > STATUS
                else
		   echo RESUBMIT, finished job $i > STATUS
		   /usr/heimdal/bin/rsh -F ekman.pdc.kth.se esubmit -n 128K -t720 /cfs/ekman/scratch/p/pschlatt/geert/run_768x385x768/runscript.sh


                endif
        else
                echo CRASHED job $i > STATUS
        endif

endif


echo Script finished.
echo
exit 0




