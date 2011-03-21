\name{FOCUS_2006_D}
\Rdversion{1.1}
\alias{FOCUS_2006_D}
\docType{data}
\title{
Dataset D from the FOCUS Kinetics report from 2006
}
\description{
Data taken from an \dQuote{existing dataset}. At each time point two data
points were generated, as the experiment was performed in duplicate.
}
\usage{data(FOCUS_2006_D)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{\code{t}}{a numeric vector containing time points}
    \item{\code{parent}}{a numeric vector containing parent concentrations 
      in percent of applied radioactivity}
    \item{\code{m1}}{a numeric vector containing concentrations of metabolite 1
    in percent of applied radioactivity}
  }
}
\source{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\examples{
data(FOCUS_2006_D)
}
\keyword{datasets}
