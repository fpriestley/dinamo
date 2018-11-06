FC = /usr/bin/gfortran
FFLAGS = -O3
LIBS = -L/home/fdp/Documents -llapack -lblas

SOURCE = main.f90 radfields.f90 readdata.f90 interpolate.f90 linsolve.f90 integrate.f90 matrixelements.f90 creategrids.f90 solve.f90 writeoutput.f90 emission.f90 particleheating.f90 particle_mod.f90 readinput.f90 countlines.f90 enthalpy.f90 calcmie.f90 bhmie.f90

main: $(SOURCE) constants_mod.o particle_mod.o
	$(FC) $(FFLAGS) -o dinamo $(SOURCE) $(LIBS)

constants_mod.o: constants_mod.f90
	$(FC) $(FFLAGS) -c constants_mod.f90

particle_mod.o: particle_mod.f90
	$(FC) $(FFLAGS) -c particle_mod.f90

.PHONY: clean

clean:
	/bin/rm *.o *.mod dinamo
