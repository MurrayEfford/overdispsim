\name{sim.cohesion}
\alias{sim.cohesion}

\title{Simulate Cohesion}

\description{

Simulates detection of a clustered population with cohesion of detection sensu
Bischof et al. (2020). The main use is to replace \code{\link{sim.capthist}} 
in the \pkg{secrdesign} function \code{\link{run.scenarios}}, as called by 
\code{\link{run_all}}.

}

\usage{

sim.cohesion(traps, popn, detectfn = 0, detectpar = list(), noccasions = 5, 
    savepopn = FALSE, gamma = 1, ...)

}

\arguments{
  \item{traps}{secr traps object}
  \item{popn}{popn object output by \code{\link[secr]{sim.popn}}}
  \item{detectfn}{integer or character value for shape of \code{\link[secr]{detectfn}}}
  \item{detectpar}{list of detection parameter values}
  \item{noccasions}{integer number of occasions}
  \item{savepopn}{logical passed to \code{\link[secr]{sim.capthist}}}
  \item{gamma}{numeric within-cluster cohesion 0-1 (Bischof et al. 2020)}
  \item{\dots}{other arguments passed to \code{\link[secr]{sim.capthist}}}
}

\details{

The argument gamma controls the degree of cohesion.

The popn argument should be a clustered population with saved attributes 
`parents' and `parentid'.

}

\value{

A capthist object.

}

\references{

Bischof, R., Dupont, P., Milleret, C., Chipperfield, J., and Royle, J. A. (2020)
Consequences of ignoring group association in spatial capture--recapture analysis.
\emph{Wildlife Biology} wlb.00649. \doi{10.2981/wlb.00649}

}

\seealso{

\code{\link{run_all}}

}
