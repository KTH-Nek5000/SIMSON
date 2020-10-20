# Example bla.i file
# Unforced spatial FSC Blasius boundary layer
#
20070716
bl0.u
bl2000_mf.u
2000.           Final time (tmax)
8500000         Max iterations (maxit)
4200000         Max CPU time (cpumax)
50000           Max wall-clock time (wallmax)
.false.         Save intermediate files (write_inter)
0.0             Time step (dt)
4               3/4 number of stages (nst)
0.8		Max CFL number (cflmaxin)
500.            Box length (xl)
.false.         Variable size (varsiz)
.0              Rotation rate (rot)
.false.         Constant mass flux (cflux)
0.              Re_tau (no effect for boundary layers)
.false.         Perturbation mode (pert)
110             Boundary condition number (ibc)
.false.         Chebyshev integration method (cim)
.false.         Galilei transformation (gall)
.false.         Constant wall-suction (suction)
.true.          Spatial simulation (spat)
.false.         Tabulated free-stream velocity (tabfree)
.false.         Read base flow from file (rbfl)
0.4             Fringe strength (fmax)
-150.           Fringe strength (fstart)
0               Fringe end (fend)
80.             Fringe rise distance (frise)
20.             Fringe fall distance (ffall)
0.000           Amplitude of oblique waves (ampob)
0.0             Amplitude of 2D-wave (amp2d)
.false.         OS-mod (osmod)
.false.         Streak (streak)
.false.         Waves (waves)
.false.         Read sgs.i (sgs)
0               SFD (isfd)
0               MHD (imhd)
0               Type of volume forcing for localized disturbance (loctyp)
.false.         Volume forcing of trip (tripf)
0               Type of wall boundary condition at lower wall (wbci)
12              cfl calc interval (icfl)
40              amp calc interval (iamp)
blkoll.amp
.false.         Generate amplitude to be written in file namamp (longli)
0               Extremum calc interval (iext)
0               xy statistics calculation interval (ixys)
.false.         Save time series (serf)
4               Number of saved 3D fields (msave)
100.            Time to save 3D field
bl100.u
300.
bl300.u
1000.
bl1000.u
1500.
bl1500.u
0               Number of wavenumbers to save (mwave)
0               Number of planes to save (npl)
