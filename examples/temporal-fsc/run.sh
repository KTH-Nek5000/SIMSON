#!/bin/sh

# Script to run the spatial Falkner-Skan-Cooke example with cross flow vortices
echo '####################################################################'
echo 'Running the temporal Falkner-Skan-Cooke example'
echo '####################################################################'

export OMP_NUM_THREADS='1'

# Generate similarity solution fsc.dat
echo 'Generating similarity solution...'
../../bls/fsc > out1.dat

# Generate initial velocity field bl0.u without disturbance
echo 'Generating undisturbed base flow...'
cp bls.noforc.i bls.i
../../bls/bls > out2.dat

# Generate initial velocity field bl0.u with disturbance
echo 'Generating base flow with eigenmode included...'
cp bls.eigenmode.i bls.i
../../bls/bls > out3.dat

# Create converged laminar flow field with eigenmode
echo 'Generating initial converged velocity field with eigenmode...'
cp bla.init.i bla.i
../../bla/bla > out4.dat

# Reset time of converged laminar flow
echo -e '1 \n bl200.u \n 1.0 \n y \n bl200_t0.u \ny\n 0.0' > cmp-parameters.dat
rm bl200_t0.u
../../cmp/cmp < cmp-parameters.dat > out5.dat
rm cmp-parameters.dat

# Run simulation with eigenmode
echo 'Running simulation...'
cp bla.forcing.i bla.i
../../bla/bla > out6.dat

echo '####################################################################'
echo 'Finished'
echo '####################################################################'
