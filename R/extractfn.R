################################################################################
# package overdisp
# extractfn. R
################################################################################

# ASSUMES environment 'local' is in workspace

extract_n <- function (ch, distribution = "poisson") {
	pop <- attr(ch, "popn")
	# Lambda <- covariates(attr(pop, 'Lambda'))$Lambda
	# if (!is.null(Lambda)) {
	# 	# ASSUMES pd applies unchanged when Lambda is a mask  UNSAFE IF DIM CHANGE
	# 	localD <- sum(.local$pd * Lambda) / sum(.local$pd)
	# }
	# else {
	localD <- NA
	tr <- traps(ch)
	# dodge fast proximity collapse to 1 occasion
	usage(tr) <- NULL
	if (distribution == 'binomial') {
		a <- sum(.local$pd) * attr(.local$mask, 'area')
		localD <- sum(pdot(pop, tr, 
						   detectpar = .local$detectpar, 
						   detectfn = 'HHN', noccasions = .local$noccasions)) / a
	}
	else {
		
		msk <- attr(pop, 'Lambda')
		if (is.null(msk)) {
			msk <- attr(pop, 'mask')
			Lambda <- covariates(msk)$D
		}
		else {
			Lambda <- covariates(msk)$Lambda
		}
		if (!is.null(Lambda)) {  # saved randomDensity
			pd <- pdot(msk, tr, detectfn = 'HHN', 
					   detectpar = .local$detectpar, 
					   noccasions = .local$noccasions)
			localD <- sum(pd * Lambda) / sum(pd)
		}
		else {
			localD <- .local$N / maskarea(msk)   # global density
		}
	}
	N <- nrow(pop)   # may be saved in run.scenarios det.args = list(savepopn=TRUE)
	if (is.null(N)) N <- NA
	n <- nrow(ch)
	chatnk <- Fletcher.chat(nk(ch), .local$enk, 3, type = 'both', verbose = FALSE)
	c(n = n, chatnk = chatnk, N = N, localD = localD)
}

# for fitted models (M)
extract_M <- function (fit, chatmin = 1) {
	if (inherits(fit, 'secr')) {
		pop <- attr(fit$capthist, "popn")
		Lambda <- covariates(attr(pop, 'Lambda'))$Lambda
		if (!is.null(Lambda)) {
			# ASSUMES pd applies unchanged when Lambda is a mask
			trueD <- sum(.local$pd * Lambda) / sum(.local$pd)
		}
		else {
			tr <- traps(fit$capthist)
			# dodge fast proximity collapse to 1 occasion
			usage(tr) <- NULL
			if (fit$details$distribution == 'binomial') {
				a <- sum(.local$pd) * attr(fit$mask, 'area')
				trueD <- sum(pdot(pop, tr, 
								  detectpar = .local$detectpar, 
								  detectfn = 'HHN', noccasions = .local$noccasions)) / a
			}
			else {
				msk <-attr(pop, 'mask')
				maskD <- covariates(msk)$D
				if (!is.null(maskD)) {  # saved randomDensity
					pd <- pdot(msk, tr, detectfn = 'HHN', 
							   detectpar = .local$detectpar, 
							   noccasions = .local$noccasions)
					trueD <- sum(pd * maskD) / sum(pd)
				}
				else {
					trueD <- .local$N / maskarea(fit$mask)
				}
			}
		}
		N <- nrow(pop)
		n <- nrow(fit$capthist)
		pred <- predict(fit)
		chatnk <- chat.nk(fit, type = 'both', verbose = FALSE)
		predF <- adjustVarD(fit, chat = chatnk['Fletcher'], chatmin = chatmin )
		predW <- adjustVarD(fit, chat = chatnk['Wedderburn'], chatmin = chatmin)
		list(N = N, n = n, chatnk = chatnk, pred = pred, predF = predF, 
			 predW = predW, trueD = trueD)
	}
	else {
		list(N = NA, n = NA, chatnk = c(NA,NA), pred = NULL, predF = NULL, 
			 predW = NULL, trueD = NA)
	}
}

extract_MCL <- function (fit) {
	if (inherits(fit, 'secr')) {
		
		pop <- attr(fit$capthist, "popn")
		Lambda <- covariates(attr(pop, 'Lambda'))$Lambda
		if (!is.null(Lambda)) {
			trueD <- sum(.local$pd * Lambda) / sum(.local$pd)
		}
		else {
			trueD <- sum(.local$pd * covariates(attr(pop, 'mask'))$D) / sum(.local$pd)
		}
		N <- nrow(pop)
		n <- nrow(fit$capthist)
		chatnk <- chat.nk(fit, type = 'both', verbose = FALSE)
		derived <- derived(fit, se.esa = TRUE)
		list(N = N, n = n, chatnk = chatnk, derived = derived, trueD = trueD)
	}
	else {
		list(N = NA, n = NA, chatnk = c(NA,NA), derived = NULL, trueD = NA)
	}
}
