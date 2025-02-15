\name{summary}
\alias{summary_n}
\alias{summary_M}
\alias{summary_MCL}

\title{Summarise Simulations}

\description{

Custom summary functions for output from \code{\link{run}}.

}

\usage{

summary_n(sims, Evar = NULL)
summary_M(sims, true = NULL)
summary_MCL(sims, true = NULL, truea = 0.0635256)

}

\arguments{

  \item{sims}{output from \code{run}}
  \item{Evar}{expected variance (should be provided for fixed-N(A))}
  \item{true}{true density}
  \item{truea}{true value of effective sampling area a}

}

\details{

To follow

}

\value{

Table with one row per scenario.

}

\seealso{

\code{\link{run}}

}


