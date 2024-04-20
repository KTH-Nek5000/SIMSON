********************************************************************************
README file for the channel flow example

Welcome to the first example in Simson!
=======================================

This examples replicates the famous case by Kim, Moin and Moser (JFM, 1987),
with a reduced resolution of 64x65x64 grid points in a domain 4pi x 2 x 2pi.

The case is run from 0 (with imposed noise) to a time of 500, with
statistics taken from 200-500.

There are examples of simple postprocessing included, which can be
run with Matlab or octave:

> pamp.m

pamp.m plots the evolution of the wall shear stress as function of time

> pstat.m

pstat.m plots the mean flow in inner and outer units.


Compiling and running instructions
==================================

1. Compile bla, bls, rit, pxyst and lambda2, by going into each directory and type
> make
Be sure to have the right resolution in the par.f file (e.g. 64x65x64 points).

2. Copy the binaries bla, bls, rit, pxyst and lambda2 into this run directory.

3. run bls to generate the initial field init.u

4. run bla to do the simulation

5. explore instantaneous velocity fields t*.u using rit

6. explore the time evolution of the wall shear stress (amp.stat) by running pamp.m in either octave or matlab

7. explore time-averaged statistics (xy.stat) by running pstat.m in either octave or matlab

8. Compute lambda2 for an instantaneous field

> lambda2 t10.u t10.u

and then run paraview to load the generated t10.u.vtr. You may use the provided channel.pvsm if needed.

That's it!

********************************************************************************