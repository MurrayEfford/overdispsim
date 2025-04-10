\name{summary}
\alias{summary_n}
\alias{summary_M}
\alias{summary_MCL}

\title{Summarise Simulations}

\description{

Custom summary functions for output from \code{\link{run_all}}.

}

\usage{

summary_n(sims, binomial = FALSE)
summary_M(sims, true = NULL)
summary_MCL(sims, true = NULL, truea = 0.0635256)

}

\arguments{

  \item{sims}{output from \code{run}}
  \item{Evar}{expected variance (should be provided for fixed-N(A))}
  \item{true}{true density}
  \item{truea}{true value of effective sampling area a}

}

\value{

Table with one row per scenario. Columns are defined as follows -

summary_n

\tabular{lll}{
[var, scale etc.] \tab Parameters specific to the scenarios \cr
varration \tab Ratio var(n) / mean(n) \cr
chatF \tab Mean 'Fletcher-chat' for count \eqn{n_k} (individuals per detector) \cr
chatW \tab Mean 'Wedderburn-chat' for count \eqn{n_k} (individuals per detector) \cr
Nsim \tab number of datasets simulated \cr
N \tab mean number of individuals in simulated population \cr
varN \tab variance of number of individuals in simulated population \cr
n \tab mean number of individuals detected \cr
varn \tab variance of number of individuals detected \cr
localD \tab mean of detection-weighted density \cr
varlocalD \tab variance of detection-weighted density \cr
VRlocalD \tab variance ratio from localD : (varlocalD/localD^2 + 1/En) / (1/En) \cr
}

summary_M

\tabular{lll}{
[var, scale etc.] \tab Parameters specific to the scenarios \cr
n \tab mean number of individuals detected \cr
N \tab mean number of individuals in simulated population \cr
nvalid \tab number of successful simulations \cr
estimate \tab mean estimated density \cr
SE.estimate \tab mean SE of estimate \cr
RSE \tab ratio of preceding \cr
trueD \tab true density; either detection-weighted (default) or global as specified by the 'true' argument\cr
RB \tab estimated relative bias relative to trueD \cr
seRB \tab SE of RB \cr
COV \tab unadjusted coverage of 95\% interval relative to trueD \cr
COVF \tab adjusted coverage \cr
chatF \tab mean 'Fletcher-chat' for count \eqn{n_k} (individuals per detector) \cr
varration \tab Ratio var(n) / mean(n) \cr
}

summary_MCL

\tabular{lll}{
[var, scale etc.] \tab Parameters specific to the scenarios \cr
n \tab mean number of individuals detected \cr
nvalid \tab number of successful simulations \cr
estimate \tab mean estimated density \cr
SE.estimate \tab mean SE of estimate \cr
RSE \tab ratio of preceding \cr
trueD \tab true density; either detection-weighted (default) or global as specified by the 'true' argument\cr
RB \tab estimated relative bias relative to trueD \cr
COV \tab unadjusted coverage of 95\% interval relative to trueD \cr
chatF \tab mean 'Fletcher-chat' for count \eqn{n_k} (individuals per detector) \cr
varration \tab Ratio var(n) / mean(n) \cr
a \tab mean a-hat effective sampling area \cr
SEa \tab SE of a-hat \cr
RBa \tab RB(a-hat) \cr
RSEa \tab RSE(a-hat) \cr
COVa \tab coverage of 95\% CI for a-hat \cr
pCVn \tab fraction of var(D-hat) attributable to var(n)\cr

}


}

\seealso{

\code{\link{run_all}}

\code{\link{plotnsum}}

}


