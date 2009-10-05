\name{kinobjects}
\Rdversion{1.1}
\alias{kinobjects}
\title{
Creates list of objects, each representing a kinetic experiment or trial
}
\description{
Function to initialise several objects representing a kinetic experiment or trial at once.
}
\usage{
  kinobjects(parent, type, systems, layers = NA, sampling_times = NA)
}
\arguments{
  \item{parent}{ The name of the parent compound }
  \item{type}{ The type of experiment or trial, optionally with an ID }
  \item{systems}{ An array of the system names }
  \item{layers}{ Optional specification of the layer names in a field trial. }
  \item{sampling_times}{ Optional specification of the sampling time pionts. }
}
\value{
  A list of lists containing the specified information.
}
\author{ Johannes Ranke }
\examples{
ko <- kinobjects("Compound XY",
	"Degradation in the environment",
	c("System 1", "System 2", "System 3"))
}
\keyword{ manip }
