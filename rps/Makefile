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
main  = rps
image = $(main)
fftp  = cvecfft_acc

# classified subroutines.
plot1 = plot1
read  = rseq.o
get   = getxys.o getxzs.o getyzs.o getxzc.o
corr  = twoxs.o twozs.o
diff  = rdcheb.o
stat  = statpl.o prob.o stagpr.o dislen.o fofx.o tmean.o intwgt.o espec.o
stat1 = accu.o
plot  = cont6.o plot1.o
raw   = raw.o wraw.o wpgmr.o frecl.o
dep   = $(main).o $(read) $(get) $(corr) $(diff) $(stat)
idep  = $(plot) $(raw) $(stat1) $(fftp).o
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
