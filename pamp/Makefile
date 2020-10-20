# ***********************************************************************
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
# ***********************************************************************

# Define fortran files in use
main    = pamp
plot    = plot1

FC = ifort
COPTS = -WB -w -r8

.f.o:
	@echo "Compiling $*.f"; $(FC) $(COPTS) -c $*.f

all: $(main)1 $(main)2 $(main)w $(main)w2

$(main)1: $(main)1.o Makefile $(plot).o
	@echo " " ; echo " Link object files into executable image '$(main)1'."
	@echo " "
	$(FC) $(COPTS) -o $@ $(main)1.o $(plot).o
	@echo " " ; echo " Done!" ; echo " "

$(main)2: $(main)2.o Makefile $(plot).o
	@echo " " ; echo " Link object files into executable image '$(main)2'."
	@echo " "
	$(FC) $(COPTS) -o $@ $(main)2.o $(plot).o
	@echo " " ; echo " Done!" ; echo " "

$(main)w: $(main)w.o Makefile $(plot).o
	@echo " " ; echo " Link object files into executable image '$(main)w'."
	@echo " "
	$(FC) $(COPTS) -o $@ $(main)w.o $(plot).o
	@echo " " ; echo " Done!" ; echo " "

$(main)w2: $(main)w2.o Makefile $(plot).o
	@echo " " ; echo " Link object files into executable image '$(main)w2'."
	@echo " "
	$(FC) $(COPTS) -o $@ $(main)w2.o $(plot).o
	@echo " " ; echo " Done!" ; echo " "

clean:
	@echo " " ; echo " Clean up all object files and all '$(main)' versions."
	@-rm -f *.o $(main)1 $(main)2 $(main)w $(main)w2 ; echo " "
