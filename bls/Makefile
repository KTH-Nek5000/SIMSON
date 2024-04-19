# ***********************************************************************
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
# ***********************************************************************

# compiler flags
#F90   = ifort
#FLAGS = -r8 -i4 -O3 -warn truncated_source -warn unused
F90   = gfortran
FLAGS = -std=legacy -fdefault-real-8 -O3


# main program and FFTs
main  = bls
image = $(main)
fftp  = cvecfft_acc

# classified subroutines
bl    = rbla.o cubip.o step.o dstep.o xlim.o
bld   = blvel.o 
noise = stnois.o
rw    = wdiscpbl.o rdiscp.o
get   = getxyp.o getxzp.o
put   = putxzp.o putxyp.o
dep   = $(main).o $(rw) $(get) $(put) $(noise)
idep  = stokes.o ran2.o icheb.o $(bl) $(bld) $(fftp).o
obj   = $(dep) $(idep)

all: fsc $(main)

# compilation rules
%.o: %.f
	 $(F90) $(FLAGS) -c $*.f 

$(image): $(obj) Makefile
	@echo " " ; echo " Link object files into executable image '$(image)'."
	@echo " "
	$(F90) $(FLAGS) -o $@ $(obj)
	@echo " " ; echo " Done!" ; echo " "

$(dep):	par.f

clean:
	@echo " " ; echo " Clean up all object files and '$(image)' and 'fsc'."
	@-rm -f *.o *.F *~ $(image) fsc ; echo " "

# Falkner-Skan-Cooke solver
fsc: fsc.o
	$(F90) $(FLAGS) -o $@ $@.o


