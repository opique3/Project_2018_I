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

#cp ../EIA_dynamics_main.f90 .
#cp ../check_MB.f90 .
#cp ../rdf.f90 .
#cp ../Makefile .
#cp ../pbc_module.f90 .
#cp ../lj_module.f90 .
#cp ../moment_module.f90 .
#cp ../position_module.f90 .
#cp ../velocities_module.f90 .
#cp ../andersen_therm_module.f90 .
#cp ../vel_verlet_module.f90 .
#cp ../print_positions_module.f90 .
#cp ../kinetic_energy_module.f90 .
#cp ../read_data_module.f90 .
#cp ../print_data_module.f90 .

make

gfortran pbc_module.f90 rdf.f90 -o rdf
gfortran check_MB.f90 -o check_MB

./dynamics input.dat
./check_MB velocity.out
./rdf traj.xyz 25

rm -f *.mod *.o dynamics check_MB rdf

