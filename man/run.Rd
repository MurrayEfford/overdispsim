\name{run}
\alias{run_all}

\title{Run Simulations}

\description{

A wrapper for \pkg{secrdesign} function \code{run.scenarios()}.

}

\usage{
run_all(nrepl, popargs, CH.function = "sim.capthist", detargs = list(savepopn = TRUE), 
    fit = FALSE, distribution = c("poisson", "binomial"), CL = FALSE, 
    start = NULL, byscenario = FALSE, extractfn    = NULL, seed = 12345, ...)
}

\arguments{
  \item{nrepl}{integer number of replicates}
  \item{popargs}{list of named arguments to
  \code{\link[secr]{sim.popn}}}
  \item{CH.function}{character name of function to simulate capthist}
  \item{detargs}{list of named arguments to
  \code{\link[secr]{sim.capthist}}}
  \item{fit}{logical; model fitted if TRUE}
  \item{distribution}{character distribution of n in fitted model ("poisson" or "binomial")}
  \item{CL}{logical; if TRUE fitting maximises the conditional likelihood}
  \item{start}{'start' argument for secr.fit}
  \item{byscenario}{logical (see secrdesign)}
  \item{extractfn}{function to extract statistics from each replicate}
  \item{seed}{integer random seed}
  \item{\dots}{additional arguments for extractfn}
}

\details{

Arguments generally follow \code{\link[secrdesign]{run.scenarios}}.

The number of threads is determined by .local$maxncores (see code).

By default, 'extractfn' is chosen depending on fit and CL:

\tabular{lll}{
fit \tab CL \tab extractfn \cr
FALSE \tab ------- \tab extract_n \cr
TRUE \tab FALSE \tab extractfn_M \cr
TRUE \tab TRUE \tab extractfn_MCL \cr
}

}

\value{

See \code{\link[secrdesign]{run.scenarios}}. Depends on extractfn.

}

\seealso{

\code{\link{extractfn}}

}
