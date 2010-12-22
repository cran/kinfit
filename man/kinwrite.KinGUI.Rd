\name{kinwrite.KinGUI}
\Rdversion{1.1}
\alias{kinwrite.KinGUI}
\encoding{latin1}
\title{
Function to write KinGUI input files from kinetic data
}
\description{
This function takes an R object as used by the \code{kinfit} package
and tries to write a text file which is compatible with the KinGUI
software tool.
}
\usage{
kinwrite.KinGUI(kinobject, file, comment=NA)
}
\arguments{
  \item{kinobject}{
	A list containing the following elements:
	The name of the parent compound to be output (\code{parent}),
	the type of the test system (\code{type}), 
	the name of the specific test system used for generating this dataset
        (\code{system}), 
	a data frame containing the raw data (\code{data}),
        which should be in the same form as required by
        by \code{\link{kinfit}}.
}
  \item{file}{
The filename, potentially including the full path, specifying where the output whould be written.
}
  \item{comment}{
An optional comment that will be integrated in the header of the KinGUI input file.
}
}
\value{
The function is called for its side effect, namely the generation of a text file.
}
\references{
Schäfer D, Mikolasch M, Rainbird P and Harvey B (2007). KinGUI: a new kinetic software tool for evaluations according to FOCUS degradation kinetics. In: Del Re AAM, Capri E, Fragoulis G and Trevisan M (Eds.). Proceedings of the XIII Symposium Pesticide Chemistry, Piacenza, 2007, p. 916-923.
}
\note{
The KinGUI software tool was announced to be freely available on CD from the main author of the paper cited above, Dieter Schäfer <dieter.schaefer@bayercropscience.com>.
}
\author{
Johannes Ranke
}

\keyword{ IO }
