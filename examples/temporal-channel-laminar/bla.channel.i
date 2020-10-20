# Example bla.i file
# Periodic laminar channel flow
#
20070716        Version of input file
ch0.u
ch40.u
40.0            Maximum simulation time (tmax)
-40             Maximum iterations (maxit)
-200            Maximum CPU time (cpumax)
-200            Maximum wall-clock time (wallmax)
.true.          Save intermediate files (write_inter)
0               Time step (dt)
4               1/3/4 number of RK stages (nst)
0.8             CFL number (cflmaxin)
48.             Box length (xl)
.false.         Variable size (varsiz)
.0              Rotation rate (rot)
.true.          Constant mass flux (cflux)
.false.         Perturbation mode (pert)
0               Boundary condition number (ibc)
.false.         Chebyshev integration method (cim)
.false.         Galilei transformation (gall)
.false.         Constant wall-suction (suction)
.false.         Spatial simulation (spat)
0.              Temporal development vel. (cdev)
.false.         Read sgs.i (sgs)
0               SFD (isfd)
0               MHD (imhd)
0               Type of localized disturbance (loctyp)
.false.         Trip forcing (tripf)
0               Boundary condition at lower wall (wbci)
4               CFL calculation interval (icfl)
4               Amplitude calculation interval (iamp)
blkoll.amp
.false.         Write urms velocities (fileurms) 
.false.         y-dependent statistics (longli)
0               Extremum calculation interval (iext)
32              xy-statistics calculation interval (ixys)
xy.stat
64              xy-statistics saving interval (ixyss)
500.            Start time for statistics (txys)
.false.         Two-point correlations (corrf)
.false.         Time series (serf)
3               Number of 3D fields to save (msave)
10.
ch10.u
20.
ch20.u
30.
ch30.u
0               Number of wavenumbers to save (mwave)
0               Number of planes to save (npl)










