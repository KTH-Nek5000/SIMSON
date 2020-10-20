init.u
450.        Re based on delta star and U freestream
600.        xl in delta star
30.         yl in delta star
34.         zl in delta star
fsc.dat
6           Flow type (Spatial Blasius) (fltype)
-80.        Base flow blending start (bstart)
80.         Base flow blending length (bslope)
.false.     Perturbation mode (pert)
0.          Shift velocity of box (ushift)
.false.     Localized disturbance (locdi)
.false.     Gaussian perturbation (gaussian)
.false.     Waves (waves)
.false.     OS modes (os)
.false.     Spectral modes (specm)
.false.     Peturbations from file (pertfromfile)
.false.     Noise (noise)
-2.         scalar 1 (m1=0 --> Dirichlet b.c. at wall)
1.                   (m1=0 --> Dirichlet b.c. at freestream)
1.          scalar 2 (m1=0.5 --> Neumann b.c. at wall)
4.                   (m1=0.5 --> Dirichlet b.c. at freestream)
-3.         scalar 3 (m1=0 --> Dirichlet b.c. at wall)
0.                   (m1=0 --> Dirichlet b.c. at freestream)
-1.         scalar 4 (m1=0.5 --> Neumann b.c.)
-2.                  (m1=0.5 --> Dirichlet b.c. at freestream)

