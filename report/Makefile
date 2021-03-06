# Generic make file for LaTeX: requires GNU make
#
# This makefile provides four targets: dvi, ps, pdf and clean.
# The default is "pdf".
# To make a dvi file, type "make dvi"
# To make a ps file, type "make ps".
# To make a pdf file, type "make pdf" or simply "make".
# To remove all files generated by make, type "make clean".
#
# Rouben Rostmaian <rostamian@umbc.edu>
# September 2002
#
# September 2005: Added the -o flag to dvips.  Thanks
# to Juan J. Merelo Guervós <jmerelo@geneura.ugr.es>
#
# Downloaded at:
# http://www.math.umbc.edu/~rouben/prosper/makefile.html

TEXFILE	:= geog5804.proj.tex

.SUFFIXES: .tex .pdf
.PHONY: all pdf clean

hostname=$(shell hostname)
ifeq ($(hostname),pontus.cee.carleton.ca)
PDF_CMD := /usr/local/texlive/2010/bin/x86_64-darwin/pdflatex --shell-escape
VIEW_CMD := open
else
PDF_CMD := pdflatex --shell-escape
VIEW_CMD := evince
endif

all: lclean ${TEXFILE:.tex=.pdf} view

${TEXFILE}: p.appendix.fixioapi_tool.tex p.eval.tex p.cover_letter.tex p.define_grid.tex p.intro.tex

pdf: ${TEXFILE:.tex=.pdf}

.tex.pdf:
	${PDF_CMD} $<;
	while grep -q "Rerun to get cross-references right." $(<:.tex=.log); \
	do \
		${PDF_CMD} $<; \
	done;

view: ${TEXFILE:.tex=.pdf}
	${VIEW_CMD} $(TEXFILE:.tex=.pdf)

clean:
	@rm -f \
	$(TEXFILE:.tex=.aux) \
	$(TEXFILE:.tex=.log) \
	$(TEXFILE:.tex=.out) \
	$(TEXFILE:.tex=.dvi) \
	$(TEXFILE:.tex=.pdf) \
	$(TEXFILE:.tex=.bib) \
	$(TEXFILE:.tex=.ps)

lclean:
	rm -f ${TEXFILE:.tex=.pdf}
