# Example bla.i file
# Spatial Blasius boundary layer with 4 passive scalars
#
20070716        Version of input file
init.u
end.uu
5000.00         Maximum simulation time (tmax)
-400             Maximum iterations (maxit)
-200            Maximum CPU time (cpumax)
-200            Maximum wall-clock time (wallmax)
.false.         Save intermediate files (write_inter)
0               Time step (dt)
4               1/3/4 number of RK stages (nst)
0.8             CFL number (cflmaxin)
600.            Box length (xl)
.false.         Variable size (varsiz)
.0              Rotation rate (rot)
.false.         Constant mass flux (cflux)
0.              Re_tau (no effect for boundary layers)
.false.         Perturbation mode (pert)
110             Boundary condition number (ibc)
0               Boundary condition for scalar 1
-2                (Dirichlet lower wall)
1.                (Dirichlet freestream)
1               Boundary condition for scalar 2
1.                (Neumann lower wall)
4.                (Dirichlet freestream)
2               Boundary condition for scalar 3
-3.               (Dirichlet lower wall)
0.                (Neumann freestream)
3               Boundary condition for scalar 4
-1.               (Neumann lower wall)
0.                (Neumann freestream)
-2.               (fixed freestream value)
.false.         Chebyshev integration method (cim)
.false.         Galilei transformation (gall)
.false.         Constant wall-suction (suction)
.true.          Spatial simulation (spat)
.false.         Tabulated free-stream velocity (tabfre)
.false.         Read base flow from file (rbfl)
1               Fringe strength (fmax)
-70.            Fringe start (fstart)
0.              Fringe end (fend)
50.             Fringe rise distance (frise)
15.             Fringe fall distance (ffall)
0.0             Amplitude oblique waves (ampob)
0.0             Amplitude 2D TS waves (amp2d)
.false.         Forcing with Orr-Sommerfeld modes (osmod)
.false.         Forcing with streaks (streak)
.false.         Forcing with waves (waves)
.false.         Read sgs.i (sgs)
0               SFD (isfd)
0               MHD (imhd)
0               Type of localized disturbance (loctyp)
# 7               Type of localized disturbance (loctyp)
# 0.1             Amplitude
# 100.            x position
# 20.             y position
# 5.              x scale
# 5.              y scale
# 0.              Frequency
.true.          Volume forcing of trip (tripf)
0.0             Steady forcing amplitude (tamps)
0.2             Time-dependent forcing amplitude (tampt)
4.0             x length scale of trip (txsc)
10.             x origin of trip (tx0)
1.0             y length scale of trip (tysc)
10              Number of z modes in trip nzt)
4.0             Time scale of trip (tdt)
-1              Random number seed for trip (seed)
0               Boundary condition at lower wall (wbci)
4               CFL calculation interval (icfl)
0               Amplitude calculation interval (iamp)
.false.         y-dependent statistics (longli)
0               Extremum calculation interval (iext)
64              xy-statistics calculation interval (ixys)
xy.stat
256             xy-statistics saving interval (ixyss)
1500.           Start time for statistics (txys)
.true.          Two-point correlations (corrf)
corr.stat
6               Number of correlations (ncorr)
.true.          Time series (serf)
ser.stat
6               Number of time series (nser)
9               Number of 3D fields to save (msave)
10.
t10.u
100.
t100.u
200.
t200.u
400.
t400.u
800.
t800.u
1600.
t1600.u
2000.
t2000.u
5000.
t5000.u
10000.
t10000.u
0          Number of wavenumbers to save (mwave)
0          Number of planes to save (npl)
endstring










