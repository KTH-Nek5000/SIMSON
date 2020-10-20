# ***********************************************************************
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
# ***********************************************************************

# usage: make [hdf=yes]

# Compiler flags
#F90     = ifort
F90     = gfortran


# add -r8 for full double precision
#FLAGS   = -O3 -fpp -openmp -warn truncated_source -warn unused
FLAGS   = -cpp -fopenmp -fdefault-real-8 -O3


# Main program and FFTs
main    = lambda2
main2   = hdf2bla
image   = $(main)
image2  = $(main2)
fftp    = cvecfft_acc
libhdf  = 
pre     =

ifeq "$(hdf)" "yes"
	pre = -DHDF
	libhdf  = -lmfhdf -ldf -ljpeg -lz
endif


# Classified subroutines
obj     = $(main).o dcheb.o eigenvalue.o comp_lam.o clenshaw.o time_string.o\
	  write_dx.o write_hdf.o write_vtk.o write_engold.o write_xml.o $(fftp).o
hdfobj  = hdf2bla.o $(fftp).o

all: $(main) $(main2)

# Compilation rules
%.o : %.f
	$(F90) $(FLAGS) $(pre) -c $*.f

$(image): $(obj) Makefile
	@echo " " ; echo " Link object files into executable image '$(image)'."
	@echo " "
	$(F90) $(FLAGS) $(pre) -o $@ $(obj) $(libhdf)
	@echo " " ; echo " Done!" ; echo " "

# Conversion to bla-format from HDF
$(image2): $(hdfobj) Makefile
	@echo " " ; echo " Link object files into executable image '$(image2)'."
	@echo " "
	$(F90) $(FLAGS) $(pre) -o $@ $(obj) $(libhdf)
	@echo " " ; echo " Done!" ; echo " "

clean:
	@echo " " ; echo " Clean up all object files, '$(image)' and '$(image2)'."
	@-rm -f *.o *.F *~ $(image) $(image2) ; echo " "

