#!/bin/sh

echo '####################################################################'
echo 'Running the spatial Falkner-Skan-Cooke example'
echo 'This will take a while...'
echo '####################################################################'

export OMP_NUM_THREADS='1'

# Generate similarity solution
echo 'Generating similary solution...'
../../bls/fsc  > out0.dat

# Generate initial velocity field bl0.u
echo 'Generating base flow...'
../../bls/bls > out1.dat

# Create converged laminar flow field without perturbations
echo 'Generating converged unforced Navier-Stokes solution...'
cp bla.init.i bla.i
../../bla/bla > out2.dat

# Change time stamp on initial velocity field
echo 'Changing time stamp on unforced Navier-Stokes solution...'
echo -e '1 \n bl2000_mf.u \n 1.0 \n y \n bl2000_mf_t0.u \ny\n 0.0' > cmp-parameters.dat
rm bl2000_mf_t0.u
../../cmp/cmp < cmp-parameters.dat > out3.dat
rm cmp-parameters.dat

# Run simulations with disturbance in the upstream part of the domain
echo 'Running simulation with trip forcing...'
cp bla.forcing.i bla.i
../../bla/bla > out4.dat

echo '####################################################################'
echo 'Finished'
echo '####################################################################'
