# Example bla.i file
# Periodic channel flow with constant mass flux
#
20070716        Version of input file
init.u
end.u
500.00        Maximum simulation time (tmax)
-40             Maximum iterations (maxit)
-200            Maximum CPU time (cpumax)
-200            Maximum wall-clock time (wallmax)
.false.          Save intermediate files (write_inter)
0               Time step (dt)
4               1/3/4 number of RK stages (nst)
0.8             CFL number (cflmaxin)
12.56637061     Box length (xl)
.true.         Variable size (varsiz)
.0              Rotation rate (rot)
.true.          Constant mass flux (cflux)
# .false.         Constant mass flux (cflux)
# 180.            Re_tau
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
16              Amplitude calculation interval (iamp)
amp.stat
.false.
.false.         y-dependent statistics (longli)
0               Extremum calculation interval (iext)
32              xy-statistics calculation interval (ixys)
xy.stat
64              xy-statistics saving interval (ixyss)
200.            Start time for statistics (txys)
.false.         Two-point correlations in z (corrf)
.false.         Two-point correlations in x (corrf_x)
.false.         Time series (serf)
5               Number of 3D fields to save (msave)
10.
t10.u
50.
t50.u
100.
t100.u
200.
t200.u
500.
t500.u
0               Number of wavenumbers to save (mwave)
0               Number of planes to save (npl)
endstring










