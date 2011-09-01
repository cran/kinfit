\name{kinfit}
\Rdversion{1.1}
\alias{kinfit}
\title{
Fit kinetic models to chemical degradation data
}
\description{
  A selection of kinetic models as defined in the FOCUS kinetics report from
  2006 are fitted to a given dataframe of chemical degradation data.
}
\usage{
kinfit(kindata, kinmodels = c("SFO"), parent.0.user = NA, parent.0.fixed = FALSE,
  start.SFO = list(parent.0 = NA, k = NA), 
  start.FOMC = list(parent.0 = NA, alpha = NA, beta = NA), 
  start.DFOP = list(parent.0 = NA, k1 = NA, k2 = NA, g = NA),
  start.HS = list(parent.0 = NA, k1 = NA, k2 = NA, tb = NA),
  algorithm = "default")
}
\arguments{
  \item{kindata}{
    A data frame containing a time variable named \code{t} and concentration 
    data for the parent compound named \code{parent}.
  }
  \item{kinmodels}{
    An array of character strings which are names of the models to be fit.
    Possible names are \code{\link{SFO}}, \code{\link{FOMC}}, \code{\link{DFOP}}
    and \code{\link{HS}}.
  }
  \item{parent.0.user}{
    The user can give a starting estimate for parent.0 here, overriding other
    potential sources for starting values as specified below.
  }
  \item{parent.0.fixed}{
    If TRUE, the initial value for the observed value will be fixed to the 
    value specified by \code{parent.0.user}.
  }
  \item{start.SFO}{
    A list of starting parameters for fitting the \code{\link{SFO}} model,
    containing \code{parent.0} and \code{k}. 
  }
  \item{start.FOMC}{
    A list of starting parameters for fitting the \code{\link{FOMC}} model,
    containing \code{parent.0}, \code{alpha} and \code{beta}. 
  }
  \item{start.DFOP}{
    A list of starting parameters for fitting the \code{\link{DFOP}} model,
    containing \code{parent.0}, \code{k1}, \code{k2} and \code{g}. 
  }
  \item{start.HS}{
    A list of starting parameters for fitting the \code{\link{HS}} model,
    containing \code{parent.0}, \code{k1}, \code{k2} and \code{tb}. 
  }
  \item{algorithm}{
    The algorithm to use for the calls to \code{\link{nls}}.
  }
}
\details{
  Per default all starting parameters are \code{NA} and the function tries to
  find suitable starting parameters on its own.
}
\value{
  A list of models of class \code{nls} representing the models that were
  fitted successfully.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\author{ Johannes Ranke }
\examples{
data(FOCUS_2006_A)
(kinfits <- kinfit(FOCUS_2006_A))
data(FOCUS_2006_B)
(kinfits <- kinfit(FOCUS_2006_B, kinmodels=c("SFO","FOMC")))
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
