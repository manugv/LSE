F90 = ifort

FFLAGS = -O3 -mcmodel=medium -shared-intel -c

LIBS =  

OBJS =  mod_precision.o mod_geometry.o mod_compact.o mod_fftw3.o mod_output.o mod_routines.o mod_lse_coeff.o main.o

corr : $(OBJS)
	$(F90) -O3 -mcmodel=medium -shared-intel -mkl $(OBJS) -o $@ $(LIBS)

mod_precision.o : mod_precision.f90
	$(F90) $(FFLAGS) mod_precision.f90

mod_geometry.o : mod_geometry.f90
	$(F90) $(FFLAGS) mod_geometry.f90

mod_compact.o :mod_compact.f90
	$(F90) $(FFLAGS) mod_compact.f90

mod_fftw3.o : mod_fftw3.f90
	$(F90) $(FFLAGS) mod_fftw3.f90

mod_output.o : mod_output.f90
	$(F90) $(FFLAGS) mod_output.f90

mod_routines.o : mod_routines.f90
	$(F90) $(FFLAGS) mod_routines.f90

mod_lse_coeff.o : mod_lse_coeff.f90
	$(F90) $(FFLAGS) mod_lse_coeff.f90

main.o : main.f90
	$(F90) $(FFLAGS) main.f90

clean:
	rm -f $(OBJS) $(OBJ) *.mod corr
