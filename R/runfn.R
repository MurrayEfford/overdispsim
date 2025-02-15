run_all <- function (
		nrepl, 
		popargs, 
		CH.function  = "sim.capthist", 
		detargs      = list(savepopn = TRUE), 
		fit          = FALSE, 
		distribution = c("poisson", "binomial"),
		CL           = FALSE,
		start        = NULL,
		byscenario   = FALSE,
		extractfn    = NULL,
		seed         = 12345,
		...
		) {
	distribution <- match.arg(distribution)
	scen <- make.scenarios (
		D          = .local$D, 
		popindex   = 1:length(popargs), 
		detectfn   = 'HHN',
		lambda0    = .local$lambda0, 
		sigma      = .local$sigma, 
		noccasions = .local$noccasions)
	
	ncores <- .local$maxncores
	# if byscenario more cores throws an error
	if (byscenario) ncores <- min(.local$ncores, nrow(scen))
	if (!fit && !byscenario) ncores <- 1   # no use for more
	
	fitargs <- list(
		detectfn = 'HHN', 
		mask     = .local$mask,
		CL       = CL, 
		start    = start,   
		details  = list(distribution = distribution))
	if (is.null(extractfn)) {
		extractfn <- if (fit) {if (CL) extract_MCL else extract_M} 
	             else extract_n
	}
	run.scenarios(
		nrepl       = nrepl, 
		scenarios   = scen, 
		trapset     = .local$traps, 
		maskset     = .local$mask, 
		pop.args    = popargs, 
		CH.function = CH.function,
		det.args    = detargs,
		fit.args    = fitargs,
		fit         = fit, 
		seed        = seed,
		extractfn   = extractfn,
		ncores      = ncores, 
		byscenario  = byscenario,
		...)
}
