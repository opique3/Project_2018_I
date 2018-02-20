opt: dynamics.f90
	gfortran -c dynamicsModule.f90
	gfortran -o dynamics dynamicsModule.o dynamics.f90 -O3

check: dynamics.f90
	gfortran -c dynamicsModule.f90
	gfortran -Wall -Wextra -fbounds-check -o dynamics dynamics.f90

backup: dynamics.f90
	rm -rf ./backup
	mkdir ./backup
	cp dynamics.f90 ./backup/
	cp dynamics ./backup/

opt-run:
	gfortran -c dynamicsModule.f90
	gfortran -o dynamics dynamicsModule.o dynamics.f90 -O3
	./dynamics input.dat

opt-mb:
	gfortran -o check_MB check_MB.f90 -O3

opt-rdf:
	gfortran -c dynamicsModule.f90
	gfortran -o rdf dynamicsModule.o rdf.f90 -O3
