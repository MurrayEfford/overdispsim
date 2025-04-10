\name{setparameters}
\alias{setparameters}

\title{Set Parameters}

\description{

Several variables are used by multiple functions and it is tidier to store them 
in an environment that is internal to the package but shared among functions. 
This function assigns values to the variables.

}

\usage{

setparameters(lambda0 = 0.5, sigma = 1, detectfn = 'HHN', noccasions = 5, 
    traps = make.grid(12, 12, detector = "proximity", spacing = 2), 
    maskspacing = 0.5, maskbuffer = 4, N = 256, maxncores = 18)
}

\arguments{
  \item{lambda0}{numeric detection parameter}
  \item{sigma}{numeric detection parameter}
  \item{detectfn}{secr \code{\link[secr]{detectfn}} code}
  \item{noccasions}{integer}
  \item{traps}{secr traps object}
  \item{maskspacing}{numeric spacing of generated mask}
  \item{maskbuffer}{numeric buffer width of generated mask}
  \item{N}{integer expected population size in masked area}
  \item{maxncores}{integer available cores for parallel processing}
}

\details{

A mask is generated. The derived variables `detectpar`, `D`, `pd` and `enk` are also computed (see vignette).

\code{setparameters} is called by `on.attach()` and may be called again to replace values.

'maxncores' is capped at the number available, as given by parallel::detectCores().

}

\value{

Invisibly returns the environment as a list with these components (order varies) -

\item{lambda0 }{intercept of HHN detection function}
\item{sigma }{spatial scale of intercept of HHN detection function}
\item{traps }{secr 'traps' object}
\item{detectpar }{list with components lambda0 and sigma*}
\item{N }{population size in mask area}
\item{D }{density*}
\item{mask }{secr 'mask' object*}
\item{maxncores }{maximum number of threads to use}
\item{noccasions }{number of sampling occasions}
\item{detectfn}{secr \code{\link[secr]{detectfn}} code}

\item{pd }{numeric overall probability of detection of AC at each mask point*}
\item{enk }{numeric vector: expected number of individuals detected at each detector*}
\item{En }{numeric expected number of individuals detected*}
\item{Evarn}{numeric expected variance of binomial n*}

* computed by setparameters(); other values are as input

}

