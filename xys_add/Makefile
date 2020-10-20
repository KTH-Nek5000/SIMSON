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
FLAGS = -fdefault-real-8 -O3

# main program, FFTs and plot routines
main  = xys_add
image = $(main)

# classified subroutines
sub = read_xys.o wxys.o
dep = $(main).o $(sub) 
obj = $(dep) 

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
