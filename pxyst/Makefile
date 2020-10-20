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
#FLAGS = -r8 -i4 -O3 -i_dynamic -mcmodel=medium -warn truncated_source -warn unused
#F90   = pgf90
#FLAGS = -r8 -O3
F90   = gfortran
FLAGS = -fdefault-real-8 -O3

# main program, FFTs and plot routines
main  = pxyst
image = $(main)
fftp  = cvecfft_acc

# classified subroutines.
plot   = plot1.o
sub    = mgrid.o rxys.o meanpl.o correlation.o correlation_old.o
comp   = pvar.o locsc.o rscale.o ddx.o ddy.o
filter = filter.o infilt.o
dep    = $(main).o $(sub) $(comp) $(filter)
idep   = cont5.o pgmr.o wpgmr.o mxpl.o mypl.o icheb.o rdcheb.o chebe2.o $(plot) $(fftp).o
obj    = $(dep) $(idep)

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
	@echo " " ; echo " Clean up all object files and '$(image)'."
	@-rm -f *.o *.F *~ $(image) ; echo " "


