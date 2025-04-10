\name{plotnsum}
\alias{plotnsum}

\title{Plot Results from summary_n}

\description{

A wrapper for commands to plot chat vs the empirical variance/mean ratio for n. 
A log scale is used for both axes.

}

\usage{

plotnsum(nsum, ny = 5, by.y = FALSE, xlim = c(0.75, 120), ylim = xlim, 
label = "", labelx = 0.15, type = c("Fletcher", "Wedderburn"), ...)

}

\arguments{

  \item{nsum}{matrix output by \code{\link{summary_n}}}
  \item{ny}{integer number of levels }
  \item{by.y}{logical; see Details}
  \item{xlim}{numeric vector of limits on x-axis}
  \item{ylim}{numeric vector of limits on y-axis}
  \item{label}{character optional label}
  \item{labelx}{numeric label offset}
  \item{type}{character type of c-hat}
  \item{\dots}{other arguments passed to points()}

}

\details{

Used for text and appendix figures in Efford & Fletcher (2024).

Each column of nsum is of length ny*nx (nx inferred). Internally the columns for 
'varration' and 'chatF' (or chatW) are formed into ny x nx matrices. Then lines 
are plotted either for each column (by.y = TRUE) or each row of these matrices.

}

\value{

none

}

\seealso{

\code{\link{summary_n}}

}

