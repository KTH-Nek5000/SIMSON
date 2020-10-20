#!/pkg/paraview/3.8.1-RC1/bin/pvpython
####/usr/local/bin/pvpython
####/usr/bin/env pvpython
#########################################################################
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
#########################################################################
#
# Driver routine for ParaView volume plot movie
#
# Note that it calls another python script from the
# command line. This is the cleanest way to do it (for now)
#
# Milos Ilak, September 16, 2010
# 
#########################################################################
import subprocess
import string
import array
import sys
from paraviewscripts import *

#
# Source data file directory
#
indir=""
gfile=""
tfile=""
ifile=""
sfiles=""

#
# Use argument as output file names, otherwise use default name
#
print "Running python driver script to generate ParaView images."
if len(sys.argv) > 1:
    print 'Arguments to driver : ' + str(sys.argv[1:])
    mfile = sys.argv[1]
else:
    print 'driver [movie_filename] (default parameter file is movie_01.dat)'
    mfile="movie_01.dat"

#
# Read global movie parameter file and generate movie*.dat file
#
indir,framefile,gfile,tfile,ifile,sfiles = setup_movie.compute_movie_globals(mfile)

#
# Create trajectory file based on the scene??.dat files
#
setup_movie.compute_orbit(sfiles,tfile,indir)
#neqlines_list = setup_movie.compute_orbit(sfiles,tfile,indir)
#sys.exit()

#
# Analyze trajectory file and find neighbouring lines based on the same dataset
#
neqlines_list = setup_movie.neqlines(tfile)

#
# Read trajectory file
#
g = open(tfile,"rt")
trajectory_data = g.readlines()
g.close()

#
# Generate grid
#
#setup_movie.gridgen(2600.0,200.0,260.0,20.0,5.0,10.0,"hairpin")
#setup_movie.gridgen(-300,5500.0,-20.0,200.0,-200.0,200.0,100.0,10.0,20.0,"hairpin")
#setup_movie.gridgen(-300,5500.0,-20.0,220.0,-200.0,200.0,200.0,40.0,40.0,"hairpin-les2")
#setup_movie.gridgen(-300,5500.0,-200.0,200.0,-20.0,200.0,100.0,20.0,10.0,"hairpin-dns")
#setup_movie.gridgen(-300,5500.0,-200.0,200.0,-20.0,220.0,200.0,40.0,40.0,"hairpin-dns2")
#sys.exit()

#
# Loop over lines in camera trajectory data file
#
i = 0
sline = 0
for line in neqlines_list:
     
    #
    # Write subset trajectory data file for each frame
    #
    g = open("tmp_traj.dat",'w')
    g.writelines(trajectory_data[sline:sline + line])
    g.close()

#    print "::: " + ifile + " ::: " + "tmp_traj.dat" + " ::: " + gfile

    subprocess.call(["pvpython",framefile,ifile,"tmp_traj.dat",gfile])
#    subprocess.call(["pvbatch","--use-offscreen-rendering",framefile,ifile,"tmp_traj.dat",gfile])

    sline = sline + line
