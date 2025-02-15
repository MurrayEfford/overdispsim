\name{myRF}
\alias{myRF}

\title{Custom Random Field}

\description{

Call RandomFields functions to generate a log-Gaussian random field.

}

\usage{

myRF(mask, parm = list(), plt = FALSE)

}

\arguments{
  \item{mask}{secr habitat mask}
  \item{parm}{list of parameters}
  \item{plt}{logical; if TRUE a plot is generated}
}

\details{

This function is used to generate a density surface for a log-Gaussian Cox Process. If fixed (a conditional LGCP) the surface is scaled to give the requested N(A). The density surface may be used in \code{\link{sim.popn}} with model2D = "IHP". (Only the unconditional process can be simulated with the spatstat function rLGCP).

The components of `parm' are 
\tabular{lll}{
Parameter \tab Default \tab Note \cr
var \tab 1 \tab in situ Gaussian variance \cr
scale \tab 1 \tab spatial scale of covariance \cr 
model \tab 'exp' \tab covariance model \cr 
N \tab 256 \tab N(A) \cr
fixed \tab TRUE \tab conditional? \cr
		
}
}

\value{

Vector of cell-specific densities.

}

\examples{

}
