# Makefile created by mkmf 19.2.0

CPPDEFS = -Duse_libMPI -Duse_netCDF -DSPMD

OTHERFLAGS = -I../MOM6-examples/build/intel/shared/repro

MK_TEMPLATE = /lustre/lysun/models/MOM6-test/MOM6-examples//src/mkmf/templates/dt2-intel.mk
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: particles
particle_driver.o: ./particle_driver.F90 particles.o particles_framework.o particles_extra.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./particle_driver.F90
particles.o: ./particles.F90 particles_framework.o particles_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./particles.F90
particles_extra.o: ./particles_extra.F90 particles_framework.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./particles_extra.F90
particles_framework.o: ./particles_framework.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./particles_framework.F90
particles_io.o: ./particles_io.F90 particles_framework.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./particles_io.F90
SRC = ./particles_extra.F90 ./particles_framework.F90 ./particles_io.F90 ./particles.F90 ./particle_driver.F90
OBJ = particles_extra.o particles_framework.o particles_io.o particles.o particle_driver.o
clean: neat
	-rm -f .particles.cppdefs $(OBJ) particles

neat:
	-rm -f $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

particles: $(OBJ)
	$(LD) $(OBJ) -o particles -L../MOM6-examples/build/intel/shared/repro -lfms $(LDFLAGS)

