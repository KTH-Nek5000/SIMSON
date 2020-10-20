# ***********************************************************************
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
# ***********************************************************************

# compiler flags
F90   = ifort
FLAGS = -r8 -i4 -O3 -warn truncated_source -warn unused

# main program, FFTs and plot routines
main  = fou
image = $(main)
fftp  = cvecfft_acc

# classified subroutines.
read = rdiscr.o
get = getxzp.o 
pg  = pgmr.o wpgmr.o
put = putxyp.o 
plot = contfou.o  splint.o curvefit.o svdcmp.o svbksb.o \
       newton.o plot1.o savelist.o foustr.o
dep = $(main).o $(read) $(get) $(put) spline.o foubranch.o
idep = $(plot) $(fftp).o $(pg)
obj = $(dep) $(idep)

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
