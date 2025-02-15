setlocalparameters <- function (
		lambda0 = 0.5,
		sigma = 1.0,
		detectfn = 'HHN',
		noccasions = 5,
		traps = make.grid(12,12, detector = 'proximity', spacing = 2.0),
		maskspacing = 0.5,
		maskbuffer = 4,
		N = 256,
		maxncores = 18)
{
	
	mask <- make.mask(traps, spacing = maskspacing, buffer = maskbuffer)
	
	.local$lambda0 <- lambda0
	.local$sigma  <- sigma
	.local$detectfn <- detectfn
	.local$detectpar <- list(lambda0 = lambda0, sigma = sigma)
	.local$traps <- traps
	.local$mask <- mask
	.local$N <- N
	.local$D <- N/maskarea(mask)
	.local$noccasions <- noccasions
	.local$maxncores <- maxncores
	
	# probability of detection at each mask point
	.local$pd <- secr::pdot(mask, traps, detectfn = detectfn, detectpar = 
					  	.local$detectpar, noccasions = noccasions)
	# expected number of individuals at each detector
	.local$enk <- secr::Enk(D = .local$D, 
					  mask, traps, detectfn = detectfn,
					  detectpar = .local$detectpar, 
					  noccasions = .local$noccasions)
	invisible(as.list(.local))
}
