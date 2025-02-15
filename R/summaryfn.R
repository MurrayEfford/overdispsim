
# summarise data simulated without fitting model

summary_n <- function (sims, Evar = NULL) {
	onescenario <- function(x) {
		n <- x[,1]
		x[n==0,] <- NA
		if (is.null(Evar)) {
			Evar <- mean(n, na.rm = TRUE)  # assume Poisson
		}
		va    <- var(n, na.rm = TRUE)
		chatF <- mean(x[,2], na.rm = TRUE)
		chatW <- if (dim(x)[2]>2) mean(x[,3], na.rm = TRUE) else NA
		N     <- if (dim(x)[2]>3) mean(x[,4], na.rm = TRUE) else NA
		varN  <- if (dim(x)[2]>3) var (x[,4], na.rm = TRUE) else NA
		Nsim      <- sum(!is.na(x[,1]))
		localD    <- if (dim(x)[2]>4) mean(x[,5], na.rm = TRUE) else NA
		varlocalD <- if (dim(x)[2]>4) var(x[,5], na.rm = TRUE) else NA
		
		data.frame(
			varration  = va/Evar, 
			chatF      = chatF, 
			chatW      = chatW, 
			Nsim       = Nsim, 
			N          = N, 
			varN       = varN, 
			n          = mean(n, na.rm = TRUE), 
			varn       = va 
			# localD     = localD, 
			# varlocalD  = varlocalD,
			# VRlocalD   = (varlocalD/localD^2 + 1/En) / (1/En),   # Davalidation
			# VRlocalDdjf = (varlocalD/localD^2 + 1/En) / (1/En)  # Davalidation
		)
	}
	parmlevels <- sapply(sims$pop.args,'[[','details')
	parmlevels <- parmlevels[rownames(parmlevels) %in% c('var','scale','mu','scale','p','A'),]
	
	out <- do.call(rbind, lapply(sims$output, onescenario))
	cbind(t(parmlevels), out)
}

# summarise fitted models

summary_M <- function (sims, true = NULL) {
	getD <- function(x, true, Dfield) {
		D <- sapply(x, function(y) unlist(y[[Dfield]]['D',c(2:5)]))
		n <- sapply(x, '[[', 'n')
		N <- sapply(x, '[[', 'N')
		chatF <- sapply(x, function(y) y$chatnk["Fletcher"]) 
		D <- data.frame(t(D))
		cbind(D, 
			  trueD  = rep(true, length.out=length(n)), 
			  RSE    = D$SE.estimate / true, 
			  RB     = (D$estimate-true) / true, 
			  COV    = D$lcl<true & D$ucl>true, 
			  n      = n, 
			  N      = N,
			  nvalid = sum(!is.na(D$SE.estimate)), 
			  chatF  = chatF)
	}
	parmlevels <- sapply(sims$pop.args,'[[','details')
	parmlevels <- parmlevels[rownames(parmlevels) %in% c('var','scale','mu','scale','p','A'),]
	outlist <- if (is.null(sims$output)) list(sims) else sims$output
	# ad hoc selection of valid simulations
	for (i in 1:length(outlist)) outlist[[i]] <- outlist[[i]][!is.na(sapply(outlist[[i]],'[[','N'))]
	
	if (is.null(true)) {
		if (!is.null(outlist[[1]][[1]]$trueD)) {
			warning("true D is detection-weighted local density")
			true <- lapply(outlist, sapply, '[[', 'trueD')
		}
		else {
			true <- N/maskarea(.local$mask)
		}
	}
	sim  <- mapply(getD, outlist, true, MoreArgs=list(Dfield='pred'), SIMPLIFY=FALSE) 
	simmean <- sapply (sim, apply, 2, mean, na.rm=T)
	simsd <- sapply (sim, apply, 2, sd, na.rm=T)
	nrepl <- sapply(sim, function(x) sum(!is.na(x$estimate)))
	covF <- mapply(getD, outlist, true, MoreArgs=list(Dfield='predF'), SIMPLIFY=FALSE) 
	out <- rbind(
		simmean, 
		seRB      = simsd['RB',]/sqrt(nrepl),
		varration = sapply (sim, function(x) var(x$n)/mean(x$n)),
		COVF      = sapply (covF, function(x) mean(x$COV, na.rm = TRUE))
	)
	keep <- c("n", "N", "nvalid", "estimate", "SE.estimate", "RSE", 
			  "trueD", "RB", "seRB", "COV", "COVF", "chatF", "varration")
	out <- t(rbind(parmlevels, out[keep,]))
	apply(out,2,as.numeric)
	
}


summary_MCL <- function (sims, true = NULL, truea = 0.0635256) {
	getD <- function(x, true) {
		getDa <- function (der) {
			if (is.null(der)) {
				list(D=rep(NA,4), a = rep(NA,4), CVnD=rep(NA,2))
			}
			else {
				D <- der['D',c(1,2,3,4)]
				a <- der['esa',c(1,2,3,4)]
				CVnD <- der['D', c('CVn','CVD')]
				list(D = D, a = a, CVnD = CVnD)
			}
		}
		est <- lapply(lapply(x, '[[','derived'), getDa)
		
		D <- do.call(rbind, lapply(est, '[[', 'D'))
		a <- do.call(rbind, lapply(est, '[[', 'a'))
		CVnD <- do.call(rbind, lapply(est, '[[', 'CVnD'))
		pCVn <- CVnD[,1]^2 / CVnD[,2]^2 
		
		n <- sapply(x, '[[', 'n')
		chatF <- sapply(x, function(y) y$chatnk["Fletcher"]) 
		cbind(D, 
			  trueD  = rep(true, length.out=length(n)), 
			  RB     = (true - D$estimate) / true, 
			  RSE    = D$SE.estimate / D$estimate,
			  COV    = D$lcl<true & D$ucl>true, 
			  n      = n, 
			  nvalid = sum(!is.na(D$SE.estimate)), 
			  a      = a$estimate,
			  SEa    = a$SE.estimate,
			  RBa    = (truea - a$estimate) / truea,
			  RSEa   = a$SE.estimate / a$estimate,
			  COVa   = a$lcl<truea & a$ucl>truea, 
			  pCVn   = pCVn,
			  chatF  = chatF)
	}
	parmlevels <- sapply(sims$pop.args,'[[','details')
	parmlevels <- parmlevels[rownames(parmlevels) %in% c('var','scale','mu','scale','p','A'),]
	outlist <- sims$output
	if (is.null(true)) {
		if (!is.null(outlist[[1]][[1]]$trueD)) {
			warning("true D is detection-weighted local density")
			true <- lapply(outlist, sapply, '[[', 'trueD')
		}
		else  {
			true <- .local$N/maskarea(.local$mask)
		}
	}
	sim  <- mapply(getD, outlist, true, SIMPLIFY=FALSE) 
	simmean <- sapply (sim, apply, 2, mean, na.rm=T)
	out <- rbind(
		simmean, 
		varration = sapply (sim, function(x) var(x$n, na.rm=TRUE)/mean(x$n, na.rm=TRUE)),
		varn  = sapply (sim, function(x) var(x$n, na.rm=TRUE)),    # var(n)
		vara  = sapply (sim, function(x) var(x$a, na.rm=TRUE)),    # var(a-hat)
		sdD   = sapply (sim, function(x) sd(x$trueD, na.rm=TRUE)),
		meanD = sapply (sim, function(x) mean(x$trueD, na.rm=TRUE))
	)
	keep <- c("n", "nvalid", "estimate", "SE.estimate", "RSE", 
			  "trueD", "RB", "COV", "chatF", "varration",
			  "a", "SEa", "RBa", "RSEa", "COVa", "pCVn")
	out <- t(rbind(parmlevels, out[keep,]))
	apply(out,2,as.numeric)
	
}
