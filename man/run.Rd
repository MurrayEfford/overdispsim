\name{run}
\alias{run_all}

\title{Run Simulations}

\description{

A wrapper for run.scenarios.

}

\usage{
run_all(nrepl, popargs, CH.function = "sim.capthist", detargs = list(savepopn = TRUE), 
    fit = FALSE, distribution = c("poisson", "binomial"), CL = FALSE, ncores = 18, 
    start = NULL, byscenario = FALSE, extractfn    = NULL, seed = 12345, ...)
}

\arguments{
  \item{nrepl}{}
  \item{popargs}{}
  \item{CH.function}{}
  \item{detargs}{}
  \item{fit}{logical}
  \item{distribution}{}
  \item{CL}{logical; if TRUE fitting maximises the conditional likelihood}
  \item{ncores}{}
  \item{start}{}
  \item{byscenario}{}
  \item{extractfn}{}
  \item{seed}{}
  \item{\dots}{additional arguments for extractfn}
}

\details{
}

\value{
}

\seealso{
}

\examples{
}
