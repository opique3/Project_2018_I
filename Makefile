
compile-main: EIA_dynamics_main.f90
	gfortran pbc_module.f90 lj_module.f90 moment_module.f90 position_module.f90 velocities_module.f90 andersen_therm_module.f90 vel_verlet_module.f90 print_positions_module.f90 kinetic_energy_module.f90 read_data_module.f90 print_data_module.f90 EIA_dynamics_main.f90 -o dynamics

compile-rdf: rdf.f90
	gfortran pbc_module.f90 rdf.f90 -o rdf

compile-pbc: pbc_module.f90
	gfortran -c pbc_module.f90

compile-lj: lj_module.f90
	gfortran -c lj_module.f90

compile-mom: moment_module.f90
	gfortran -c moment_module.f90

compile-init_pos: position_module.f90
	gfortran -c position_module.f90

compile-init_vel: velocities_module.f90
	gfortran -c velocities_module.f90

compile-thermo: andersen_therm_module.f90
	gfortran -c andersen_therm_module.f90

compile-verlet: vel_verlet_module.f90
	gfortran -c vel_verlet_module.f90

compile-pos: print_positions_module.f90
	gfortran -c print_positions_module.f90

compile-KE: kinetic_energy_module.f90
	gfortran -c kinetic_energy_module.f90

compile-read_data: read_data_module.f90
	gfortran -c read_data_module.f90

compile-print_data: print_data_module.f90
	gfortran -c print_data_module.f90

compile-MB: check_MB.f90
	gfortran check_MB.f90 -o check_MB

