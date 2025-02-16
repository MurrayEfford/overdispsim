\name{extend}
\alias{extend}

\title{Extend Argument List}

\description{

Expand a single list of arguments using all combinations of values.

}

\usage{

extend(baseargs, values)

}

\arguments{
  \item{baseargs}{list of arguments for \code{sim.popn}}
  \item{values}{list of value vectors}
}

\details{

The vectors in 'values' refer to 'details' arguments of \code{sim.popn}.

}

\value{

List of argument lists for \code{\link{sim.popn}}.

}

\seealso{

\code{\link{setlocalparameters}}

}

\examples{

setlocalparameters()
basepopargs3 <- list(D = randomDensity, 
					 core = .local$mask,
					 buffer = 0,
					 model2D = "IHP",
					 details = list(D = .local$D, p=0.5, A=0.25, rescale=TRUE))

# vary A, p arguments of randomDensity					 
popargs3 <- extend(basepopargs3, values = list(A = c(0.25,0.5), p = c(0.25, 0.5)) )

# view details
sapply(popargs3, '[[', 'details')

}