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
FLAGS = -O3 -fdefault-real-8
# main program and FFTs
main  = cmp
image = $(main)
fftp  = cvecfft_acc

# classified subroutines
read  = rdisca.o
write = wdisca.o
get   = getxy.o getxz.o
put   = putxy.o putxz.o
amp   = wcamp.o namp.o boxamp.o
ext   = wcext.o boxext.o plnmax.o plnmin.o fft2db.o fft2df.o
comp  = add.o comxz.o dcheb.o power.o
symm  = symm_z.o
dep   = $(main).o $(read) $(write) $(get) $(put) $(comp) $(symm) $(amp) $(ext)
idep  = $(fftp).o
obj   = $(dep) $(idep)

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
