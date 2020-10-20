#!/bin/sh

echo '####################################################################'
echo 'Running the temporal Couette turbulent spot flow case'
echo 'This will take a while...'
echo '####################################################################'

export OMP_NUM_THREADS='1'

# Generate initial velocity field c0.u
echo 'Generating initial condition flow field...'
cp bls.localized.i bls.i
../../bls/bls > out1.dat

# Run simulation
echo 'Running simulation...'
cp bla.spot.i bla.i
../../bla/bla > out2.dat

echo '####################################################################'
echo 'Finished'
echo '####################################################################'
