# see also setparameters()
# environment is populated by onAttach()
.local <- new.env(parent = emptyenv())

resetchatmin <- function (sims, chatmin) {
	outlist <- if (is.null(sims$output)) list(sims) else sims$output
	onerepl <- function (out) {
		out$predF <- adjustVarD(out$pred, chatmin = chatmin, chat = out$chatnk['Fletcher'])
		out$predW <- adjustVarD(out$pred, chatmin = chatmin, chat = out$chatnk['Wedderburn'])
		out
	}	
	sims$output <-lapply(outlist, function(x) lapply(x, onerepl))
	sims
}

Blcl <- function (p, n, alpha=0.05) 
	## Binomial confidence interval
	## Wilson score interval with continuity correction after Newcombe (1998)
	## Newcombe, Robert G. "Two-Sided Confidence Intervals for the Single Proportion:
	##    Comparison of Seven Methods," Statistics in Medicine, 17, 857-872 (1998).
{
	z <- qnorm(1-alpha/2)
	L <- (2*n*p + z^2 - 1 - z * (z^2 - 2 -1/n + 4 * p *(n*(1-p) + 1))^0.5) / (2 * (n + z^2))
	ifelse (p == 0, 0, L)
}

Bucl <- function (p, n, alpha=0.05) 
	## Binomial confidence interval
	## Wilson score interval with continuity correction after Newcombe (1998)
	## Newcombe, Robert G. "Two-Sided Confidence Intervals for the Single Proportion:
	##    Comparison of Seven Methods," Statistics in Medicine, 17, 857-872 (1998).
{
	z <- qnorm(1-alpha/2)
	U <- (2*n*p + z^2 + 1 + z * (z^2 + 2 -1/n + 4 * p *(n*(1-p) - 1))^0.5) / (2 * (n + z^2))
	ifelse (p == 1, 1, U)
}

addCL <- function (x, summ, cov = 'COV', nrepl = 1000, ...) {
	# nrepl <- summ['notNA','nvalid']
	nrepl <- summ[,'nvalid']
	if (cov == 'RB') {
		segments (x, summ[,'RB']-2*summ[,'seRB'], x, summ[,'RB']+2*summ[,'seRB'])
	}
	else {
		segments (x, Blcl(summ[,cov],nrepl, ...), x, Bucl(summ[,cov],nrepl, ...))
	}
}
