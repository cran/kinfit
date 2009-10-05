R.exe -e "Sweave('kinfit.Rnw', stylepath=FALSE)"
pdflatex.exe kinfit
bibtex.exe kinfit
pdflatex.exe kinfit
pdflatex.exe kinfit
