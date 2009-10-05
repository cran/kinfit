\name{kinresplot}
\Rdversion{1.1}
\alias{kinresplot}
\title{
Creates a plot of the residual for specified kinetic fits
}
\description{
Function to create a residual plot for a specified fitted model
}
\usage{
kinresplot(kinobject, kinmodel, xlab = "Time [days]", ylab = "Residual [\% of applied radioactivity]", maxabs = "auto")
}
\arguments{
  \item{kinobject}{
	A list containing the following elements:
	The name of the parent compound to be output (\code{parent}),
	the type of the test system (\code{type}), 
	the name of the specific test system used for generating this dataset
        (\code{system}), 
	the list of fitted kinetic models (\code{fits}), as returned by
        \code{\link{kinfit}}, and optionally the list of results
        (\code{results}) as returned by \code{\link{kinresults}}.
        Also optional is the label position of the test compound (\code{label})
        and the source of the data (\code{source}). }
  \item{kinmodel}{ The fitted model for which the residuals should be plotted. }
  \item{xlab}{ Label for the x axis. }
  \item{ylab}{ Label for the y axis. }
  \item{maxabs}{ Maximum value of the absolute residuals, will be calculated
        from the residuals if not specified otherwise. }
}
\value{
The function is called for its side effect, namely creating a residual plot
for the specified fit.
}
\author{ Johannes Ranke }
\examples{
data(FOCUS_2006_C)
kinfits <- kinfit(FOCUS_2006_C)
kinobject <- list(
	parent = "Compound XY",
	type = "Degradation in the environment",
	system = "System 1",	
	source = "Synthetic example data from FOCUS kinetics",
	data = FOCUS_2006_C,
	fits = kinfits,
	results = kinresults(kinfits))	
\dontrun{kinresplot(kinobject, "SFO")}
}
\keyword{ hplot }
