# Example bla.i file
# Periodic laminar channel flow
#
20070716        Version of input file
init.u
end.uu
400.0           Maximum simulation time (tmax)
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
.false.          Constant mass flux (cflux)
0
.false.         Perturbation mode (pert)
0               Boundary condition number (ibc)
.false.         Chebyshev integration method (cim)
.false.         Galilei transformation (gall)
.false.         Constant wall-suction (suction)
.true.          Spatial simulation (spat)
.false.         Tabulated free-stream velocity (tabfre)
.true.          Read base flow from file (rbfl)
baseflow.u
1.              Fringe strength (fmax)
-13.            Fringe start (fstart)
0.              Fringe end (fend)
7.              Fringe rise distance (frise)
3.              Fringe fall distance (ffall)
0.0             Amplitude oblique waves (ampob)
0.0             Amplitude 2D TS waves (amp2d)
.false.         Forcing with Orr-Sommerfeld modes (osmod)
.false.         Forcing with streaks (streak)
.false.         Forcing with waves (waves)
.false.         Read sgs.i (sgs)
0               SFD (isfd)
0               MHD (imhd)
0               Type of localized disturbance (loctyp)
.false.         Trip forcing (tripf)
0               Boundary condition at lower wall (wbci)
4               CFL calculation interval (icfl)
0               Amplitude calculation interval (iamp)
.false.         y-dependent statistics (longli)
0               Extremum calculation interval (iext)
0               xy-statistics calculation interval (ixys)
6               Number of 3D fields to save (msave)
10.
t10.u
20.
t20.u
30.
t30.u
40.
t40.u
50.
t50.u
60.
t60.u
0               Number of wavenumbers to save (mwave)
0               Number of planes to save (npl)










