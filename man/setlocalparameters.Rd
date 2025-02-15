\name{setlocalparameters}
\alias{setlocalparameters}

\title{Set Parameters}

\description{

Several variables are used by multiple functions and it is tidier to store them in an environment that is internal to the package but shared among functions. This function assigns values to the variables.

}

\usage{

setlocalparameters(lambda0 = 0.5, sigma = 1, noccasions = 5, 
    traps = make.grid(12, 12, detector = "proximity", spacing = 2), 
    maskspacing = 0.5, maskbuffer = 4, N = 256, maxncores = 18)
}

\arguments{
  \item{lambda0}{numeric detection parameter}
  \item{sigma}{numeric detection parameter}
  \item{noccasions}{integer}
  \item{traps}{secr traps object}
  \item{maskspacing}{numeric spacing of generated mask}
  \item{maskbuffer}{numeric buffer width of generated mask}
  \item{N}{integer expected population size in masked area}
  \item{maxncores}{integer available cores for parallel processing}
}

\details{

A mask is generated. The derived variables `detectpar`, `D`, `pd` and `enk` are also computed (see vignette).

\code{setlocalparameters} is called by `on.attach()` and may be called again to replace values.

}

\value{

Invisibly returns the environment as a list with components -
\item{lambda0 }{}
\item{sigma }{}
\item{detectpar }{}
\item{noccasions }{}
\item{N }{}
\item{D }{}
\item{traps }{}
\item{mask }{}
\item{pd }{}
\item{enk }{}
\item{maxncores }{}

}

