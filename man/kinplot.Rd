\name{kinplot}
\Rdversion{1.1}
\alias{kinplot}
\title{
Creates a plot of the kinetic fits
}
\description{
Function to create a plot for a set of fitted models
}
\usage{
kinplot(kinobject, xlab = "Time [days]", ylab = "Parent [\% of applied radioactivity]", ylim = c("auto", "auto"), lpos = "topright")
}
\arguments{
  \item{kinobject}{
	A list containing the following elements:
	The name of the parent compound to be output (\code{parent}),
	the type of the test system (\code{type}), 
	the name of the specific test system used for generating this dataset
        (\code{system}), 
	the list of fitted kinetic models (\code{fits}), as returned by
        \code{\link{kinfit}}, and the list of results (\code{results}) 
        as returned by \code{\link{kinresults}}.
	Optionally also
	the label position of the test compound (\code{label}) and
	the source of the data (\code{source}). }
  \item{xlab}{ Label for the x axis. }
  \item{ylab}{ Label for the y axis. }
  \item{ylim}{ An array of length two holding the range for values on the y axis or "auto". }
  \item{lpos}{ Where should the legend be placed? Will be passed on to
		\code{\link{legend}}. }
}
\value{
The function is called for its side effect, namely creating a plot with 
the fitted model.
}
\author{ Johannes Ranke }
\examples{
data(FOCUS_2006_C)
kinfits <- kinfit(FOCUS_2006_C, kinmodels = c("SFO", "FOMC", "DFOP"))
kinobject <- list(
	parent = "Compound XY",
	type = "Degradation in the environment",
	system = "System 1",	
	source = "Synthetic example data from FOCUS kinetics",
	data = FOCUS_2006_C,
	fits = kinfits,
	results = kinresults(kinfits))	
\dontrun{kinplot(kinobject)}
}
\keyword{ hplot }
