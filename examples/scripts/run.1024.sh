#!/bin/tcsh
module add beta-modules
module add pgi/10.3
module add openmpi/1.4.2-pgi
#module add i-compilers
#module add mvapich/1.0

cd /cfs/ekman/scratch/p/pschlatt/geert/run_768x385x768

cp -f $SP_HOSTFILE old_host_file.$SP_JID

bash -c 'cat old_host_file.$SP_JID | while read line; do for ((i=1;i<=8;i+=1)); do echo "${line}"; done; done > new_host_file.$SP_JID'

mpirun --mca mpi_leave_pinned 0 --mca mpi_paffinity_alone 1 -n 1024  -machinefile ./new_host_file.$SP_JID ./blaq.1024 > bla.out.$SP_JID

