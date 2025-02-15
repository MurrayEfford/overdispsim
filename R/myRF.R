myRF <- function (mask, parm = list(), plt = FALSE) 
{
	if (!(requireNamespace("RandomFields"))) {
		stop("install RandomFields")
	}
	defaultparm <- list(
		var = 1, 
		scale = 1, 
		model = 'exp', 
		N = 256,
		fixed = TRUE)
	parm <- replace(defaultparm, names(parm), parm)
	D <- parm$N / maskarea(mask)
	mu <- log(D) - parm$var/2
	
	model <- match.arg(parm$model[1], c('exp','gauss'))
	# default is negative exponential autocorrelation
	if (model == 'exp')
		model <- RandomFields::RMexp(var = parm$var, scale = parm$scale)   
	else if (model == 'gauss')
		model <- RandomFields::RMgauss(var = parm$var, scale = parm$scale) 
	else
		stop("unknown model")
	
	model <- model + RandomFields::RMtrend(mean = mu)
	
	rf <- RandomFields::RFsimulate(model, x = as.matrix(mask))@data[,1]
	if (parm$fixed) {
		pi <- exp(rf) / sum(exp(rf))
		D  <- pi * D
	}
	else {
		D  <- exp(rf)
	}
	
	if (plt) {
		covariates(mask) <- data.frame(D = D)
		plot(mask, cov = 'D', dots = FALSE, border=2)
		parm1 <- parm[c('var','scale')]
		mtext(side = 3, line = 0, cex=0.7, paste(names(parm1), parm1, collapse='  '))
	}
	return(D)   # vector of cell-specific densities
}
