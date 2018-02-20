#!/bin/csh
#$ -N 'MD_test'
#$ -o job.out
#$ -e job.err
#$ -q cerqt2.q
#$ -S /bin/csh
#$ -cwd
#$ -pe smp 1


setenv old `pwd`
cd $old
cp ../develop/dynamicsModule.f90 .
cp ../develop/dynamics.f90 .
cp ../develop/check_MB.f90 .
cp ../develop/rdf.f90 .

gfortran -c dynamicsModule.f90
gfortran -o dynamics dynamicsModule.o dynamics.f90 -O3
gfortran -o check_MB check_MB.f90 -O3
gfortran -o rdf dynamicsModule.o rdf.f90 -O3

./dynamics input.dat
./check_MB velocity.out
./rdf traj.xyz 25

rm -f dynamics.* dynamicsModule.* dynamicsmodule* dynamics

