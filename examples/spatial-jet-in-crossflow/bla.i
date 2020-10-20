20070716
init.u
out.u
0.5	        	tmax
-200	       	maxiterations
248400     		cpumax  
248400     		wallmax
.false.    		save intermediate files
-1.0        		time step
4           		1/3/4 number of time integration stages
0.9		cflmaxin CFL number
75.        	boxlength xl
.false.     		varsiz variable size
0.          		rotation rate
.true.			const mass flux
.false.		pert Perturb. eqs.  .true.
110         		boundary condition number
.false.     		Chebyshev integration method (cim)
.false.     		galilei transformation
.false.     		suction activated
.true.      		spatial simulation (spat)
.false.     		tabfre tabulated free-stream velocity
.false.     		rbfl .true. = 3D base flow from file "nambfl" ref8.u
1.         		fmax fringe strength
-15.	         	fstart
0.           		fend
8.          		frise
3.          		ffall
0.0         		ampob ampl. of oblique waves
0.0         		amp2d ampl. of 2d ts waves
.false.      		forcing with Orr-Sommerfeld modes osamp Sqb289.fsm	osfil
.false.		Streaks
.false.		Waves
.false.		LES
0		isfd
0		mhd
0          		loctyp type of loc. volume forcing 
.false.     		volume forcing of trip (tripf) 
-2        		wbci: type of wall boundary condition at the wall  
1			jet profile (1=parabolic, 2=top hat)
3.0			jet diameter (in units of dstar)
9.375			xjet (streamwise coord., in outer units)
0.5			zjet (spanwise coord., fraction of domain width)
1.0			velocity ratio R
4           		icfl cfl-number calc interval
0           		amp calc interval (iamp)
.false.     		y-dependent statistics .false.     soloper bla_ctrl
0           		extremum calc interval (iext)
0           		xy statistics calculation interval (ixys)
0		number of saved interm. 3-d fields (msave)
0           		number of saved wavenumbers (mwave)
0           		number of saved planes (npl)
endstring
























