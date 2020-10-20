# ***********************************************************************
#
# $HeadURL$
# $LastChangedDate$
# $LastChangedBy$
# $LastChangedRevision$
#
# ***********************************************************************

# Makefile for documentation

# Gobal configuration file
include ../config.mk

# Gobal Makefile rules
include ../rules.mk

TEXFILES=\
simson-user-guide.tex \
simson-svn-guide.tex \
simson-paraview-guide.tex

PDFFILES=\
simson-user-guide.pdf \
simson-svn-guide.pdf \
simson-paraview-guide.pdf

# Printed version

#TEXFILESv4=\
#simson-user-guide-v4.0.tex

#PDFFILESv4=\
#simson-user-guide-v4.0.pdf



###########################################################################

#all: $(PDFFILES) $(PDFFILESv4)

all: $(PDFFILES)

.SUFFIXES: .pdf .tex

%.pdf: %.tex simson.cls simson.bst fluids.bib Makefile
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<
	pdflatex $<

clean:
	@echo " " ; echo " Cleaning up all auxiliary files"
	@-rm -f *.dvi *.log *.aux *.pdf *.out *.toc *.lof *.bbl *.blg ; echo " "
