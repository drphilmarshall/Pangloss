

all: sulireport

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out *.bbl paper-lcr.ps paper-lcr.pdf DELIVERABLE-SULI-REPORT-EVERETT.pdf


# ----------------------------------------------------------------------------

sulireport:
	pdflatex suli-report-swe.tex
	pdflatex suli-report-swe.tex
#	bibtex suli-report-swe
# 	pdflatex suli-report-swe.tex
#	pdflatex suli-report-swe.tex
	-mv suli-report-swe.pdf DELIVERABLE-SULI-REPORT-EVERETT.pdf

# ----------------------------------------------------------------------------

paper: paper-lcr.pdf
a4: paper-lcr-a4.pdf

# paper-lcr.pdf: paper-lcr.tex macros.tex
# 	pdflatex paper-lcr.tex
# 	pdflatex paper-lcr.tex
# 	bibtex paper-lcr
# 	pdflatex paper-lcr.tex
# 	pdflatex paper-lcr.tex

# paper-lcr.bbl: paper-lcr.tex references.bib
# 	bubble -f paper-lcr.tex references.bib >& paper-lcr.bbl.log
# 	cat paper-lcr.bbl.log

paper-lcr.dvi: paper-lcr.tex macros.tex
	latex paper-lcr.tex
	latex paper-lcr.tex
# 	bibtex paper-lcr
# 	latex paper-lcr.tex
# 	latex paper-lcr.tex

paper-lcr.pdf: paper-lcr.dvi
# 	\dvipdf -sPAPERSIZE=a4 paper-lcr.dvi
	\dvipdf -sPAPERSIZE=letter paper-lcr.dvi

paper-lcr-a4.pdf: paper-lcr.dvi
	\dvipdf -sPAPERSIZE=a4 paper-lcr.dvi paper-lcr-a4.pdf

# paper-lcr.ps: paper-lcr.dvi
# # 	\dvips -t a4 -o paper-lcr.ps paper-lcr.dvi
# 	\dvips -t letter -o paper-lcr.ps paper-lcr.dvi

# ----------------------------------------------------------------------------
