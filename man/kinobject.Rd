\name{kinobject}
\Rdversion{1.1}
\alias{kinobject}
\title{
Creates list representing a kinetic experiment or trial
}
\description{
Function to initialise an object representing a kinetic experiment or trial.
}
\usage{
  kinobject(parent, type, system, layers = NA, sampling_times = NA)
}
\arguments{
  \item{parent}{ The name of the parent compound }
  \item{type}{ The type of experiment or trial, optionally with an ID }
  \item{system}{ System name. Important if several systems were used, e.g. several soils. }
  \item{layers}{ Optional specification of the layer names in a field trial. }
  \item{sampling_times}{ Optional specification of the sampling time pionts. }
}
\value{
  A list containing the specified information.
}
\author{ Johannes Ranke }
\examples{
ko <- kinobject("Compound XY",
	"Degradation in the environment",
	"System 1")
}
\keyword{ manip }
