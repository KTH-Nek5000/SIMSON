********************************************************************************
README file for the Simson code (on github)


Welcome to Simson!
==================
Simson is a code containing a pseudo-spectral solver for the
incompressible Navier-Stokes equations together with pre- and
postprocessing tools. Both direct numerical simulations (DNS) and
large-eddy simulations (LES) are supported.


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


Root directory files of Simson
==============================
config.mk - Includes parameters read by all Makefiles
configure - Script that configures the system and stores the
            configuration in config.mk
COPYRIGHT - Copyright and contact information document
par.f     - Compile time parameters regarding resolution etc.
README    - This document
rules.mk  - Generic Makefile build rules
todo.txt  - A text file including things to correct/fix/add etc.


Installing
==========
The makefiles named Makefile.config should not be edited. Instead, the
script configure should be run. The configure script updates the
config.mk file in the root directory based on analyzes of the system
or on information from a configuration file. Note that there exists a
default config.mk file which can be directly edited and used which
means that it is not necessary to run the configure script.

  ./configure --help


Compiling
=========
Before starting to compile the code the resolution and some other
compile-time parameters have to be chosen for the flow problem that
you want to run. These parameters are fetched from the par.f file
located in the same directory as the application that is about to be
compiled. Note however that there is a template par.f in the root
directory which can be distributed to the relevant directories by
writing

  make dist

Now you can run

  make

which will compile all the programs that are part of the Simson
package. To install the files in the directory based on the $MACHINE
variable type

  make install

To install it to, for example, test/bin directory simply write

  make prefix=test install

From the main directory it is also possible to compile only a specific
tool by writing

  make rit.all

To install rit type

  make rit.install

Note that if the code is not compiled the install directive will also
do that.

Note also that in each directory there are two Makefiles, one named
Makefile.config and one named Makefile. It is the former Makefile that
is being used by the configure script setup. The other Makefile can be
edited and run locally without dependencies from config.mk etc. This
could be convenient, for example, when testing different compiler
options etc.


Running
=======
For more details on how to install, compile and run the code read the manual
located in the doc directory. Also view the examples directory.


********************************************************************************
