\name{sim.cohesion}
\alias{sim.cohesion}

\title{Simulate Cohesion}

\description{

Simulates detection of a clustered population with cohesion of detection sensu Bischof et al. (2020). The main use is to replace \code{\link{sim.capthist}} in the \pkg{secrdesign} function \code{\link{run.scenarios}}, as called by \code{\link{run}}.

}

\usage{

sim.cohesion(traps, popn, detectfn = 0, detectpar = list(), noccasions = 5, 
    savepopn = FALSE, gamma = 1, ...)

}

\arguments{
  \item{traps}{}
  \item{popn}{}
  \item{detectfn}{}
  \item{detectpar}{}
  \item{noccasions}{}
  \item{savepopn}{}
  \item{gamma}{}
  \item{\dots}{}
}

\details{

The argument gamma controls the degree of cohesion.

The popn argument should be a clustered population with saved attributes `parents' and `parentid'.

}

\value{

A capthist object.

}

\references{

Bischof, R., Dupont, P., Milleret, C., Chipperfield, J., and Royle, J. A. (2020) Consequences of ignoring group association in spatial capture--recapture analysis. *Wildlife Biology* wlb.00649. DOI 10.2981/wlb.00649

}

\seealso{

\code{\link{run}}

}

\examples{
}
