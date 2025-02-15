\name{extractfn}
\alias{extractfn}
\alias{extract_n}
\alias{extract_M}
\alias{extract_MCL}

\title{Extract Simulation Result}

\description{

Functions to extract summary statistics from each replicate.

}

\usage{

extract_n(ch, distribution = "poisson")
extract_M(fit, chatmin = 1)
extract_MCL(fit)

}

\arguments{
  \item{ch}{secr capthist object}
  \item{distribution}{character}
  \item{fit}{fitted secr model}
  \item{chatmin}{numeric threshold of c-hat for adjustment of variance}
}

\details{

The value of 'distribution' may be set to "binomial" for fixed-N(A).

}

\value{

For extract_n, a numeric vector

   \tabular{ll}{
   n \tab number of detected individuals \cr
   chatnk \tab Fletcher's c-hat \cr
   N \tab number of simulated individuals \cr
   localD \tab \cr
}

For extract_M, a list with components: 
   
   \tabular{ll}{
   N \tab number of simulated individuals \cr
   n \tab number of detected individuals \cr
   chatnk \tab Fletcher's c-hat \cr
   pred \tab table of estimates from predict() \cr
   predF \tab table of estimates with variance adustment using Fletcher's c-hat \cr
   predW \tab table of estimates with variance adustment using Wedderburn c-hat \cr
   trueD \tab detection-adjusted true density \cr
}

For extract_MCL, a list with components: 
   
   \tabular{ll}{
   N \tab number of simulated individuals \cr
   n \tab number of detected individuals \cr
   chatnk \tab Fletcher's c-hat \cr
   derived \tab table of estimates from derived() \cr
   trueD \tab detection-adjusted true density \cr
}

}

\seealso{

\code{\link{run_all}}

}

