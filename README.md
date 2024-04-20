********************************************************************************
README file for the Simson code (on public github)


Welcome to Simson!
==================
Simson is a code containing a pseudo-spectral solver for the
incompressible Navier-Stokes equations together with selected pre- and
postprocessing tools. Both direct numerical simulations (DNS) and
large-eddy simulations (LES) are supported. This version does only
include one-dimensional parallelisation.


Copyright and contact
=====================
See COPYRIGHT file.


Content of main directories of Simson
=====================================
bla       - Main program
bls       - Program to generate initial velocity (with or without heat) fields
cmp       - Program to subtract and compare velocity fields
config    - Contains configuration files for different systems
common    - Common subroutines between programs
doc       - Documentation
examples  - Example flow cases
fou       - Program to Fourier transform velocity fields in time
matlab    - Matlab scripts
pamp      - Programs to plot amplitude data from one or more amplitude files
paraview  - Scripts to generate ParaView based animations
pext1     - Program to plot a number of components from an extremum file
pxyst     - Program to do xy statistics plotting
rit       - Program to plot solutions from complete velocity fields
rps       - Program to plot planes
xys_add   - Program to compute xyz statistics
lambda2   - Program co compute lambda2, Q for visualisation


Root directory files of Simson
==============================
COPYRIGHT - Copyright and contact information document
README    - This document
todo.txt  - A text file including things to correct/fix/add etc.


Installing
==========
Compile using make in each subdirectory the various codes. Make sure to
adapt par.f in each subdirectory. For bla there are additional compilation
options such as omp, mpi and arch.

As discussed, before starting to compile the code the resolution and some other
compile-time parameters have to be chosen for the flow problem that
you want to run. These parameters are fetched from the par.f file
located in the same directory as the application that is about to be
compiled. 


Running
=======
For more details on how to install, compile and run the code read the manual
located in the doc directory. Also view the examples directory.

Usually, the stacksize of the shell needs to be extended, e.g. using

> ulimit -s unlimited

The visualisation codes typically require an xterm. Thus run rit, pxyst, etc.
from an xterm (rather than the built in terminal).


********************************************************************************
